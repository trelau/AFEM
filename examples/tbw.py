from __future__ import print_function

from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.misc.check import pairwise
from afem.structure import *
from afem.topology import *

Settings.log_to_console(True)

# Import model
fname = r'..\models\TBW_SUGAR.stp'
ImportVSP.step_file(fname)

# Get bodies.
fuselage = ImportVSP.get_body('Fuselage')
wing = ImportVSP.get_body('Wing')
other_wing = ImportVSP.get_body('Wing.2')
gear = ImportVSP.get_body('Gear Pod')
other_gear = ImportVSP.get_body('Gear Pod.1')
htail = ImportVSP.get_body('Htail')
other_htail = ImportVSP.get_body('Htail.3')
vtail = ImportVSP.get_body('Vtail')
jury = ImportVSP.get_body('Jury')
other_jury = ImportVSP.get_body('Jury.4')
strut = ImportVSP.get_body('Strut')
other_strut = ImportVSP.get_body('Strut.5')

for body in [fuselage, wing, other_wing, gear, other_gear, htail, other_htail,
             vtail]:
    body.set_color(0.5, 0.5, 0.5)
    body.set_transparency(0.5)
    Viewer.add(body)

# WING ------------------------------------------------------------------------
AssemblyAPI.create_assy('wing assy')

# Center wing structure will be based on the intersection between the wings
# and fuselage.
joined_wings = FuseShapes(wing, other_wing).shape
shape = IntersectShapes(fuselage, joined_wings).shape

# Use the y-dimensions of the bounding box of the intersection curve.
bbox = BBox()
bbox.add_shape(shape)
ymax = bbox.ymax

# Center wing box
xz_plane = PlaneByAxes().plane
root_rib_pln = PlaneByAxes((0., ymax, 0.), 'xz').plane

# Front center spar
p1 = wing.eval(0.15, 0.)
pln = PlaneByAxes(p1, 'yz').plane
fc_spar = SparBetweenShapes('fc spar', xz_plane, root_rib_pln, wing, pln).spar

# Rear center spar
p1 = wing.eval(0.70, 0.)
pln = PlaneByAxes(p1, 'yz').plane
rc_spar = SparBetweenShapes('rc spar', xz_plane, root_rib_pln, wing, pln).spar

# Root rib
root_rib = RibByPoints('root rib', fc_spar.p2, rc_spar.p2, wing).rib

# Strut rib where strut intersects.
p0 = strut.eval(0.5, 1.)
pln = PlaneByAxes(p0, 'xz').plane
strut_rib = RibBySurface('kink rib', pln, wing).rib

# Tip rib
tip_rib = RibByParameters('tip rib', 0.15, 0.99, 0.70, 0.99, wing).rib

# Inboard front spar
u = strut_rib.local_to_global_u(0.15)
p2 = strut_rib.point_on_cref(u)
pln = PlaneByIntersectingShapes(root_rib, fc_spar, p2).plane
inbd_fspar = SparByPoints('inbd fspar', root_rib.p1, p2, wing, pln).spar

# Inboard rear spar
u = strut_rib.local_to_global_u(0.70)
p2 = strut_rib.point_on_cref(u)
pln = PlaneByIntersectingShapes(root_rib, rc_spar, p2).plane
inbd_rspar = SparByPoints('inbd rspar', root_rib.p2, p2, wing, pln).spar

# Outboard front spar
pln = PlaneByIntersectingShapes(strut_rib, inbd_fspar, tip_rib.p1).plane
outbd_fspar = SparByPoints('outbd fspar', inbd_fspar.p2, tip_rib.p1, wing,
                           pln).spar

# Outboard rear spar
pln = PlaneByIntersectingShapes(strut_rib, inbd_rspar, tip_rib.p2).plane
outbd_rspar = SparByPoints('outbd rspar', inbd_rspar.p2, tip_rib.p2, wing,
                           pln).spar

# Jury rib where strut intersects
p0 = jury.eval(0.5, 1.)
pln = PlaneByAxes(p0, 'xz').plane
jury_rib = RibBetweenShapes('jury rib', inbd_fspar, inbd_rspar, wing, pln).rib

# Inboard ribs
u2 = inbd_rspar.invert_cref(jury_rib.p2)
inbd_ribs = RibsAlongCurveByDistance('inbd rib', inbd_rspar.cref, 30.,
                                     inbd_fspar, inbd_rspar, wing, u2=u2,
                                     d1=18, d2=-30).ribs

# Middle ribs
mid_ribs = RibsAlongCurveByDistance('mid rib', inbd_rspar.cref, 30.,
                                    inbd_fspar, inbd_rspar, wing, u1=u2,
                                    d1=18, d2=-30).ribs

# Outboard ribs
outbd_ribs = RibsAlongCurveByDistance('outbd rib', outbd_rspar.cref, 30.,
                                      outbd_fspar, outbd_rspar, wing,
                                      d1=18, d2=-30).ribs

# Join and discard internal structure.
internal_parts = AssemblyAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

# Wing skin
wskin = SkinByBody('wing skin', wing).skin

# Join the wing skin and internal structure
wskin.fuse(*internal_parts)

# Discard faces touching reference surface.
wskin.discard_by_dmin(wing.sref_shape, 1.)

# Fix skin since it's not a single shell anymore, but a compound of two
# shells (upper and lower skin).
wskin.fix()

Viewer.add(*AssemblyAPI.get_parts())

# HTAIL -----------------------------------------------------------------------
AssemblyAPI.create_assy('htail assy')

fspar = SparByParameters('htail fspar', 0.15, 0.01, 0.15, 0.99, htail).spar
rspar = SparByParameters('htail rspar', 0.70, 0.01, 0.70, 0.99, htail).spar
root = RibByPoints('htail root rib', fspar.p1, rspar.p1, htail).rib
tip = RibByPoints('htail tip rib', fspar.p2, rspar.p2, htail).rib
ribs = RibsAlongCurveByDistance('htail rib', rspar.cref, 12, fspar, rspar,
                                htail, d1=12, d2=-12).ribs

internal_parts = AssemblyAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

skin = SkinByBody('htail skin', htail).skin
skin.fuse(*internal_parts)
skin.discard_by_dmin(htail.sref_shape, 1.)
skin.fix()

Viewer.add(*AssemblyAPI.get_parts())

# VTAIL -----------------------------------------------------------------------
AssemblyAPI.create_assy('vtail assy')

u, v = ProjectPointToSurface(fspar.p1, vtail.sref).nearest_param
mspar = SparByParameters('vtail mspar', u, 0.05, u, v, vtail).spar

u, v = ProjectPointToSurface(rspar.p1, vtail.sref).nearest_param
rspar = SparByParameters('vtail rspar', u, 0.05, u, v, vtail).spar
fspar = SparByParameters('vtail fspar', 0.12, 0.05, 0.12, v, vtail).spar
RibByPoints('vtail root rib', fspar.p1, rspar.p1, vtail)
RibByPoints('vtail tip rib', fspar.p2, rspar.p2, vtail)
RibsAlongCurveByDistance('vtail rib', rspar.cref, 18., fspar, rspar, vtail,
                         d1=18, d2=-72)

internal_parts = AssemblyAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

skin = SkinByBody('vtail skin', vtail).skin
skin.fuse(*internal_parts)
skin.discard_by_dmin(vtail.sref_shape, 1.)
skin.fix()

Viewer.add(*AssemblyAPI.get_parts())

# FUSELAGE --------------------------------------------------------------------
AssemblyAPI.create_assy('fuselage assy')

pln = PlaneByAxes((60, 0, 0), 'yz').plane
bh1 = BulkheadBySurface('bh 1', pln, fuselage).bulkhead

pln = PlaneByAxes((180, 0, 0), 'yz').plane
bh2 = BulkheadBySurface('bh 2', pln, fuselage).bulkhead
p0 = fc_spar.p1
pln = PlaneByAxes(p0, 'yz').plane
bh3 = BulkheadBySurface('bh 3', pln, fuselage).bulkhead

p0 = rc_spar.p1
pln = PlaneByAxes(p0, 'yz').plane
bh4 = BulkheadBySurface('bh 4', pln, fuselage).bulkhead

p0 = fspar.p1
pln = PlaneByAxes(p0, 'yz').plane
bh5 = BulkheadBySurface('bh 5', pln, fuselage).bulkhead

p0 = mspar.p1
pln = PlaneByAxes(p0, 'yz').plane
bh6 = BulkheadBySurface('bh 6', pln, fuselage).bulkhead

p0 = rspar.p1
pln = PlaneByAxes(p0, 'yz').plane
bh7 = BulkheadBySurface('bh 7', pln, fuselage).bulkhead

pln = PlaneByAxes((0, 0, -24), 'xy').plane
floor = FloorBySurface('floor', pln, fuselage).floor

shell = ShellByFaces(bh1.faces).shell
fwd = HalfspaceByShape(shell, (-1e6, 0, 0)).solid
floor.cut(fwd)
shell = ShellByFaces(bh5.faces).shell
aft = HalfspaceByShape(shell, (1e6, 0, 0)).solid
floor.cut(aft)

plns = [bh1.sref, bh2.sref, bh3.sref, bh4.sref, bh5.sref, bh6.sref, bh7.sref]

indx = 1
frames = []
for pln1, pln2 in pairwise(plns):
    builder = FramesBetweenPlanesByDistance('frame', pln1, pln2, 24., fuselage,
                                            3.5, first_index=indx)
    frames += builder.frames
    indx = builder.next_index

internal_parts = AssemblyAPI.get_parts()

skin = SkinByBody('fuselage skin', fuselage).skin
skin.set_transparency(0.5)

FuseSurfaceParts([skin], internal_parts)

Viewer.add(*AssemblyAPI.get_parts())

# MIRROR ----------------------------------------------------------------------
xz_pln = PlaneByAxes().plane

parts = AssemblyAPI.get_parts('wing assy')
for part in parts:
    part.set_mirror(xz_pln)

parts = AssemblyAPI.get_parts('htail assy')
for part in parts:
    part.set_mirror(xz_pln)

# GEAR ------------------------------------------------------------------------
AssemblyAPI.create_assy('gear assy')

p0 = gear.eval(0.15, 1.)
pln = PlaneByAxes(p0, 'yz').plane
spar1 = SparBySurface('gear spar 1', pln, gear).spar

p0 = gear.eval(0.70, 1.)
pln = PlaneByAxes(p0, 'yz').plane
spar2 = SparBySurface('gear spar 2', pln, gear).spar

plns = PlanesAlongCurveByDistance(spar2.cref, 30, d1=30, d2=-30).planes
i = 1
for pln in plns:
    name = ' '.join(['gear rib', str(i)])
    RibBySurface(name, pln, gear)
    i += 1

internal_parts = AssemblyAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)

skin = SkinByBody('gear skin', gear).skin
skin.fuse(*internal_parts)
skin.set_transparency(0.5)

# STRUTS ----------------------------------------------------------------------
crv = strut.extract_curve(0.5, 0., 0.5, 1.)
crv.set_color(1, 0, 0)
beam1 = BeamByCurve('beam 1', crv).beam
beam1.set_color(1, 0, 0)

crv = jury.extract_curve(0.5, 0., 0.5, 1.)
crv.set_color(1, 0, 0)
beam2 = BeamByCurve('beam 2', crv).beam
beam2.set_color(1, 0, 0)

p2 = beam1.p1
ProjectPointToSurface(p2, wing.sref, update=True)
beam3 = BeamByPoints('beam 3', beam1.p1, p2).beam
beam3.set_color(1, 0, 0)

# VIEW ------------------------------------------------------------------------
parts = AssemblyAPI.get_parts()
for part in parts:
    part.set_mirror(xz_pln)

Viewer.add(*AssemblyAPI.get_parts())
Viewer.show()
