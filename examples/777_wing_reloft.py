from afem.config import Settings
from afem.exchange import ImportVSP
from afem.geometry import *
from afem.graphics import Viewer
from afem.topology import *

v = Viewer()

Settings.log_to_console()

# Do not split or divide closed surfaces
fn = r'../models/777-200LR_nosplit.stp'
ImportVSP.step_file(fn, divide_closed=False)
wing = ImportVSP.get_body('Wing')

wing.set_transparency(0.5)

# Reference curve along sref
cref = wing.sref.u_iso(0.25)

# Points where cuts will be along cref. Offset by 1 at root and tip because
# the xz-plane may not perfectly intersect the wing sections since they may be
# angled differently. Could prepend/append isocurves at root and tip?
cref = NurbsCurve.downcast(cref)
builder = PointsAlongCurveByDistance(cref, 60., d1=1., d2=-1.)
pnts = builder.points
parms = builder.parameters

# Uncomment this to include knot locations of sref
# knots = cref.knots
# parms += list(knots)
# parms.sort()

# At each point find the intersection of the wing and xz-plane. For some
# reason, OCC is struggling when using BOP intersection. Instead,
# use a common operation to find the face, then extract the outer wire. This
# will be reported as a bug.
wires = []
for u in parms:
    p = cref.eval(u)
    xz_pln = PlaneByAxes(p, 'xz').plane
    xz_face = FaceBySurface(xz_pln).face
    bop = CommonShapes(wing.solid, xz_face)

    # Build wire(s) from the shapes
    builder = WiresByShape(bop.shape)
    assert builder.nwires == 1
    wire = builder.wires[0]
    assert wire.Closed()
    v.add(wire)
    wires.append(wire)

v.add(*wires)
v.start()
v.view.remove_all()

# Loft the sections into a continuous surface. This may result in unexpected
# behavior since you are potentially skinning a complex shape with a single
# surface.
loft = LoftShape(wires, False, check_compatibility=True)

# Extract the face from the loft. Hopefully there is only one.
faces = ExploreShape.get_faces(loft.shape)
assert len(faces) == 1
face = faces[0]

# Extract the surface from the face
srf = ExploreShape.surface_of_face(face)

# Should be a NURBS surface
srf = NurbsSurface.downcast(srf)

v.add(srf)
v.start()
v.view.remove_all()

srf.set_color(0.5, 0.5, 0.5)
srf.set_transparency(0.5)

# Extract an iso curve at v-parameter in u-direction
c = srf.v_iso(0.5)
c = NurbsCurve.downcast(c)

v.add(srf)
v.add(c)
v.start()
v.view.remove_all()

# Grab a column of control points and move them to see what happens
i1, i2 = srf.locate_v(0.5, with_knot_repetition=True)
print(i1, i2)

cp = srf.cp
for i in range(1, srf.n + 1):
    # Use "- 1" since OCC index starts at 1 while Python starts at 0
    p = Point(*cp[i - 1, i1 - 1])
    v.add(p)
    new_p = p.copy()
    new_p.translate((0., 0., 48.))
    v.add(new_p)
    srf.set_cp(i, i1, new_p)
v.add(srf)
v.start()
v.view.remove_all()

# Open issues:
# - Figure out right index for bounding control points
# - How to smoothly move all control points

# Can access raw VSP surface this way and modify control points directly
vsp_surf = wing.get_metadata('vsp surface')

v.add(vsp_surf)
v.start()
v.view.remove_all()

# Grab some control points and modify
cp = vsp_surf.cp
indx = 3
for i in range(1, vsp_surf.m + 1):
    p = Point(*cp[indx - 1, i - 1])
    v.add(p)
    new_p = p.copy()
    new_p.translate((0., 0., 48.))
    v.add(new_p)
    vsp_surf.set_cp(indx, i, new_p)
v.add(vsp_surf)
v.start()
