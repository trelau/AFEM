from asap.config import Settings
from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP, StepExport
from asap.structure import AssemblyData, CreatePart, PartTools
from asap.topology import ShapeTools

# Set units to inch.
Settings.set_units('in')

# Import model
fn = './models/GT_supersonic_delta_split_nomerge.stp'
ImportVSP.step_file(fn)

# Get components.
wing = ImportVSP.get_body('WING')
fuselage = ImportVSP.get_body('FUSELAGE')

skin = CreatePart.skin.from_body('skin', wing)
print(ShapeTools.get_tolerance(skin))

bodies = ImportVSP.get_bodies()
for name in bodies:
    body = bodies[name]
    body.set_transparency(0.5)
    body.set_color(0.5, 0.5, 0.5)
    Viewer.add_items(body)

# OUTBOARD WING ---------------------------------------------------------------
# Get spanwise parameter of reference surface at section 2. If that small
# section is NOT removed, then this should be wing.vknots[2].
v = wing.vknots[1]

# Rib at section 2
outbd_rib = CreatePart.rib.by_parameters('outbd root rib', wing, 0.12, v,
                                         0.80, v)

# Tip rib oriented parallel to outboard edge of wing. Use derivative of
# reference surface to define a plane.
p0 = wing.sref.eval(0.5, 0.99)
vnorm = wing.sref.deriv(0., 1., 0, 1)
pln = CreateGeom.plane_by_normal(p0, vnorm)
outbd_tip_rib = CreatePart.rib.by_sref('outbd tip rib', wing, pln)

# Outbd front spar
p1 = outbd_rib.p1
u = outbd_tip_rib.local_to_global(0.15)
p2 = outbd_tip_rib.eval(u)
outbd_fspar = CreatePart.spar.by_points('outbd front spar', wing,
                                        p1, p2)

# Outbd rear spar
p1 = outbd_rib.p2
u = outbd_tip_rib.local_to_global(0.80)
p2 = outbd_tip_rib.eval(u)
outbd_rspar = CreatePart.spar.by_points('outbd rear spar', wing,
                                        p1, p2)

# Outbd mid spar
u = outbd_rib.local_to_global(0.5)
p1 = outbd_rib.eval(u)
u = outbd_tip_rib.local_to_global(0.5)
p2 = outbd_tip_rib.eval(u)
outbd_mspar = CreatePart.spar.by_points('outbd mid spar', wing,
                                        p1, p2)

# Ribs along outbd root rib fwd and aft of mid spar.
p0 = outbd_rspar.p1
vnorm = outbd_rspar.cref.deriv(0., 1)
ref_pln = CreateGeom.plane_by_normal(p0, vnorm)
u2 = outbd_rib.invert(outbd_mspar.p1)
ribs = CreatePart.rib.along_curve('outbd rib', wing, outbd_rib.cref,
                                  outbd_fspar.sref, outbd_rib.sref, 24.,
                                  s1=24., s2=-12., ref_pln=ref_pln, u2=u2)

indx = len(ribs) + 1
u2 = outbd_rib.invert(outbd_mspar.p1)
ribs += CreatePart.rib.along_curve('outbd rib', wing, outbd_rib.cref,
                                   outbd_fspar.sref, outbd_rib.sref, 24.,
                                   s1=12, s2=-12., ref_pln=ref_pln, u1=u2,
                                   indx=indx)

# Outbd ribs
indx = len(ribs) + 1
CreatePart.rib.along_curve('outbd rib', wing, outbd_rspar.cref,
                           outbd_fspar.sref, outbd_rspar.sref, 24.,
                           s1=12., s2=-24., indx=indx)

# Center rib
pln = CreateGeom.plane_by_axes(axes='xz')
center_rib = CreatePart.rib.by_parameters('center rib', wing, 0., 0.,
                                          1., 0., pln)

# Front inbd spar
u = center_rib.local_to_global(0.05)
p1 = center_rib.eval(u)
p2 = outbd_rib.p1
pln = ShapeTools.plane_from_section(outbd_rib, outbd_fspar, p1)
inbd_fspar = CreatePart.spar.by_points('inbd front spar', wing, p1, p2, pln)

# Inbd rear spar
p1 = outbd_rib.p2
center_rib.project_points([p1])
p2 = outbd_rib.p2
pln = ShapeTools.plane_from_section(outbd_rib, outbd_rspar, p1)
inbd_rspar = CreatePart.spar.by_points('inbd rear spar', wing, p1, p2, pln)

# Inbd mid spar at outbd front spar intersection.
p1 = outbd_rib.p1
center_rib.project_points([p1])
p2 = outbd_rib.p2
pln = ShapeTools.plane_from_section(outbd_rib, outbd_fspar, p1)
inbd_mspar1 = CreatePart.spar.by_points('inbd mid spar 1', wing, p1, p2, pln)

# Inbd mid spar at outbd mid spar intersection.
p1 = outbd_mspar.p1
center_rib.project_points([p1])
p2 = outbd_mspar.p1
pln = ShapeTools.plane_from_section(outbd_rib, outbd_mspar, p1)
inbd_mspar2 = CreatePart.spar.by_points('inbd mid spar 2', wing, p1, p2, pln)

# Inbd forward spars
u1 = center_rib.invert(inbd_fspar.p1)
u2 = center_rib.invert(inbd_mspar1.p1)
spars = CreatePart.spar.along_curve('inbd fwd spar', wing, center_rib.cref,
                                    center_rib.sref, inbd_fspar.sref, 72,
                                    s1=72, s2=-72, u1=u1, u2=u2)

# Inbd ribs are only placed at intersection of inbd fwd spars and inbd front
# spar.
i = 1
for spar in spars:
    p1 = spar.p2
    p2 = spar.p2
    inbd_rspar.project_points([p2])
    pln = ShapeTools.plane_from_section(inbd_fspar, spar, p2)
    name = ' '.join(['inbd rib', str(i)])
    CreatePart.spar.by_points(name, wing, p1, p2, pln)
    i += 1

# Fuse internal structure
wing_parts = AssemblyData.get_parts()
PartTools.fuse_wing_parts(wing_parts)
PartTools.discard_faces(wing_parts)

Viewer.add_items(*AssemblyData.get_parts())
Viewer.show()

# Create skin. For now this replaces the center rib with wing skin since
# they overlap. Not ideal...but fix later or define center rib differently.
# skin = CreatePart.skin.from_body('skin', wing)
print('valid before fuse', ShapeTools.is_valid(skin),
      ShapeTools.get_tolerance(skin, 1))
skin.fuse(*wing_parts)
skin.set_transparency(0.5)
print(skin.fix())
print(ShapeTools.is_valid(skin))

for part in AssemblyData.get_parts():
    print(part.label, ShapeTools.get_tolerance(part))

Viewer.add_items(*AssemblyData.get_parts())
Viewer.show()

step = StepExport()
shape = ShapeTools.make_compound(AssemblyData.get_parts())
step.transfer(shape)
step.write('gt_supersonic_delta.stp')

# Try IGES export.
from OCC.IGESControl import IGESControl_Writer, IGESControl_Controller_Init


IGESControl_Controller_Init()
iges = IGESControl_Writer('IN', 1)
for part in AssemblyData.get_parts():
    print(part.label, part, part.ShapeType())
    try:
        print(iges.AddShape(part))
    except RuntimeError:
        print('Failed to translate shape to IGES file.')
print(iges.Write('gt_supersonic_dela.igs'))

# Try STL export.
from asap.io import StlExport
from OCC.BRepMesh import BRepMesh_IncrementalMesh

# Increase tessellation.
BRepMesh_IncrementalMesh(shape, 1.0e-7)

stl = StlExport()
stl.write(shape, 'gt_supersonic_delta.stl')
