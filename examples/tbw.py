from __future__ import print_function

import time

from OCC.BRepBndLib import brepbndlib_Add
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Bnd import Bnd_Box
from OCC.gce import gce_MakeMirror

from asap.fem import MeshData
from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import AssemblyData, CreatePart, PartTools
from asap.topology import ShapeTools

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

for name in ImportVSP.get_bodies():
    body = ImportVSP.get_body(name)
    body.set_color(0.5, 0.5, 0.5)
    body.set_transparency(0.5)
    Viewer.add_items(body)

all_parts = []

# WING ------------------------------------------------------------------------
AssemblyData.create_assy('wing assy')

# Center wing structure will be based on the intersection between the wings
# and fuselage.
joined_wings = wing.fuse(other_wing)
shape = ShapeTools.bsection(fuselage, joined_wings)

# Use the y-dimensions of the bounding box of the intersection curve.
bbox = Bnd_Box()
brepbndlib_Add(shape, bbox)
corner_max = bbox.CornerMax()
ymax = corner_max.Y()

# Center wing box
xz_plane = CreateGeom.plane_by_axes(axes='xz')
root_rib_pln = CreateGeom.plane_by_axes([0, ymax, 0], 'xz')

# Front center spar
p1 = wing.eval(0.15, 0.)
pln = CreateGeom.plane_by_axes(p1, 'yz')
fc_spar = CreatePart.spar.between_geom('fc spar', wing, xz_plane,
                                       root_rib_pln, pln)
# Rear center spar
p1 = wing.eval(0.70, 0.)
pln = CreateGeom.plane_by_axes(p1, 'yz')
rc_spar = CreatePart.spar.between_geom('fc spar', wing, xz_plane,
                                       root_rib_pln, pln)

# Root rib
root_rib = CreatePart.rib.by_points('root rib', wing, fc_spar.p2, rc_spar.p2)

# Strut rib where strut intersects.
p0 = strut.eval(0.5, 1.)
pln = CreateGeom.plane_by_axes(p0, 'xz')
strut_rib = CreatePart.rib.by_sref('kink rib', wing, pln)

# Tip rib
tip_rib = CreatePart.rib.by_parameters('tip rib', wing, 0.15, 0.99, 0.70, 0.99)

# Inboard front spar
u = strut_rib.local_to_global_u(0.15)
p2 = strut_rib.eval_cref(u)
pln = ShapeTools.plane_from_section(root_rib, fc_spar, p2)
inbd_fspar = CreatePart.spar.by_points('inbd front spar', wing, root_rib.p1,
                                       p2, pln)

# Inboard rear spar
u = strut_rib.local_to_global_u(0.70)
p2 = strut_rib.eval_cref(u)
pln = ShapeTools.plane_from_section(root_rib, rc_spar, p2)
inbd_rspar = CreatePart.spar.by_points('inbd rear spar', wing, root_rib.p2,
                                       p2, pln)

# Outboard front spar
pln = ShapeTools.plane_from_section(strut_rib, inbd_fspar, tip_rib.p1)
outbd_fspar = CreatePart.spar.by_points('outbd front spar', wing,
                                        inbd_fspar.p2, tip_rib.p1, pln)

# Outboard rear spar
pln = ShapeTools.plane_from_section(strut_rib, inbd_rspar, tip_rib.p2)
outbd_rspar = CreatePart.spar.by_points('outbd rear spar', wing,
                                        inbd_rspar.p2, tip_rib.p2, pln)

# Jury rib where strut intersects.
p0 = jury.eval(0.5, 1.)
pln = CreateGeom.plane_by_axes(p0, 'xz')
jury_rib = CreatePart.rib.between_geom('jury rib', wing, inbd_fspar.sref,
                                       inbd_rspar.sref, pln)

# Inboard ribs.
u2 = inbd_rspar.invert_cref(jury_rib.p2)
inbd_ribs = CreatePart.rib.along_curve('inbd rib', wing, inbd_rspar.cref,
                                       inbd_fspar.sref, inbd_rspar.sref,
                                       30., u2=u2, s1=18., s2=-30.)

# Middle ribs
mid_ribs = CreatePart.rib.along_curve('mid rib', wing, inbd_rspar.cref,
                                      inbd_fspar.sref, inbd_rspar.sref,
                                      30., u1=u2, s1=18., s2=-30.)

# Outboard ribs.
outbd_ribs = CreatePart.rib.along_curve('outbd rib', wing, outbd_rspar.cref,
                                        outbd_fspar.sref, outbd_rspar.sref,
                                        30., s1=18., s2=-30.)

# Join and discard internal structure.
internal_parts = AssemblyData.get_parts()
PartTools.fuse_wing_parts(internal_parts)
PartTools.discard_faces(internal_parts)

# Wing skin
wskin = CreatePart.skin.from_body('wing skin', wing)

# Join the wing skin and internal structure
wskin.fuse(*internal_parts)

# Discard faces touching reference surface.
wskin.discard(wing.sref)

# Fix skin since it's not a single shell anymore, but a compound of two
# shells (upper and lower skin).
wskin.fix()

# Viewer.add_items(*AssemblyData.get_parts())
all_parts += AssemblyData.get_parts()

# HTAIL -----------------------------------------------------------------------
AssemblyData.create_assy('htail assy')

fspar = CreatePart.spar.by_parameters('htail fspar', htail, 0.15, 0.01,
                                      0.15, 0.99)

rspar = CreatePart.spar.by_parameters('htail rspar', htail, 0.70, 0.01,
                                      0.70, 0.99)

root = CreatePart.rib.by_points('htail root rib', htail, fspar.p1, rspar.p1)
tip = CreatePart.rib.by_points('htail tip rib', htail, fspar.p2, rspar.p2)

ribs = CreatePart.rib.along_curve('htail rib', htail, rspar.cref,
                                  fspar.sref, rspar.sref,
                                  12, s1=12., s2=-12)

internal_parts = AssemblyData.get_parts()
PartTools.fuse_wing_parts(internal_parts)
PartTools.discard_faces(internal_parts)

skin = CreatePart.skin.from_body('htail skin', htail)
skin.fuse(*internal_parts)
skin.discard(htail.sref)
skin.fix()

# Viewer.add_items(*AssemblyData.get_parts())
all_parts += AssemblyData.get_parts()

# VTAIL -----------------------------------------------------------------------
AssemblyData.create_assy('vtail assy')

u, v = vtail.invert_point(fspar.p1)
mspar = CreatePart.spar.by_parameters('vtail mspar', vtail, u, 0.05,
                                      u, v)

u, v = vtail.invert_point(rspar.p1)
rspar = CreatePart.spar.by_parameters('vtail rspar', vtail, u, 0.05,
                                      u, v)

fspar = CreatePart.spar.by_parameters('vtail fspar', vtail, 0.075, 0.05,
                                      0.075, v)

root = CreatePart.rib.by_points('vtail root rib', vtail, fspar.p1, rspar.p1)
tip = CreatePart.rib.by_points('vtail tip rib', vtail, fspar.p2, rspar.p2)

ribs = CreatePart.rib.along_curve('vtail rib', vtail, rspar.cref,
                                  fspar.sref, rspar.sref,
                                  18, s1=18., s2=-18)

internal_parts = AssemblyData.get_parts()
PartTools.fuse_wing_parts(internal_parts)
PartTools.discard_faces(internal_parts)

skin = CreatePart.skin.from_body('vtail skin', vtail)
skin.fuse(*internal_parts)
skin.discard(vtail.sref)
skin.fix()

# Viewer.add_items(*AssemblyData.get_parts())
all_parts += AssemblyData.get_parts()

# FUSELAGE --------------------------------------------------------------------
AssemblyData.create_assy('fuselage assy')

pln = CreateGeom.plane_by_axes([60, 0, 0], 'yz')
bh1 = CreatePart.bulkhead.by_sref('bh 1', fuselage, pln)

pln = CreateGeom.plane_by_axes([180, 0, 0], 'yz')
bh2 = CreatePart.bulkhead.by_sref('bh 2', fuselage, pln)

p0 = fc_spar.p1
pln = CreateGeom.plane_by_axes(p0, 'yz')
bh3 = CreatePart.bulkhead.by_sref('bh 3', fuselage, pln)

p0 = rc_spar.p1
pln = CreateGeom.plane_by_axes(p0, 'yz')
bh4 = CreatePart.bulkhead.by_sref('bh 4', fuselage, pln)

p0 = fspar.p1
pln = CreateGeom.plane_by_axes(p0, 'yz')
bh5 = CreatePart.bulkhead.by_sref('bh 5', fuselage, pln)

p0 = mspar.p1
pln = CreateGeom.plane_by_axes(p0, 'yz')
bh6 = CreatePart.bulkhead.by_sref('bh 6', fuselage, pln)

p0 = rspar.p1
pln = CreateGeom.plane_by_axes(p0, 'yz')
bh7 = CreatePart.bulkhead.by_sref('bh 7', fuselage, pln)

pln = CreateGeom.plane_by_axes([0, 0, -24], 'xy')
floor = CreatePart.floor.by_sref('floor', fuselage, pln)

# above_floor = ShapeTools.make_halfspace(floor, [0, 0, 1e6])
# bh3.cut(above_floor)
# bh4.cut(above_floor)

fwd = ShapeTools.make_halfspace(bh1, [-1e6, 0, 0])
floor.cut(fwd)
aft = ShapeTools.make_halfspace(bh5, [1e6, 0, 0])
floor.cut(aft)

skin = CreatePart.skin.from_body('fuselage skin', fuselage)

plns = [bh1.sref, bh2.sref, bh3.sref, bh4.sref, bh5.sref, bh6.sref, bh7.sref]
frames = CreatePart.frame.between_planes('frame', fuselage, plns, 3.5, 24.)

PartTools.fuse_parts(AssemblyData.get_parts())

# Viewer.add_items(*AssemblyData.get_parts())
all_parts += AssemblyData.get_parts()

# MIRROR ----------------------------------------------------------------------


xz_pln = CreateGeom.plane_by_axes(axes='xz')

trsf = gce_MakeMirror(xz_plane.Pln()).Value()

transform = BRepBuilderAPI_Transform(trsf)

parts = AssemblyData.get_parts('wing assy')
compound = ShapeTools.make_compound(parts)
transform.Perform(compound, True)
shape = transform.Shape()
# Viewer.display(shape)
all_parts.append(shape)

parts = AssemblyData.get_parts('htail assy')
compound = ShapeTools.make_compound(parts)
transform.Perform(compound, True)
shape = transform.Shape()
# Viewer.display(shape)
all_parts.append(shape)

# GEAR ------------------------------------------------------------------------
AssemblyData.create_assy('gear assy')

p0 = gear.eval(0.15, 1.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
spar1 = CreatePart.spar.by_sref('gear spar 1', gear, pln)

p0 = gear.eval(0.70, 1.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
spar2 = CreatePart.spar.by_sref('gear spar 2', gear, pln)

plns = CreateGeom.planes_along_curve(spar2.cref, 30., s1=30., s2=-30.)
i = 1
for pln in plns:
    name = ' '.join(['gear rib', str(i)])
    CreatePart.rib.by_sref('gear rib', gear, pln)
    i += 1

internal_parts = AssemblyData.get_parts()
PartTools.fuse_wing_parts(internal_parts)
skin = CreatePart.skin.from_body('gear skin', gear)
skin.fuse(*internal_parts)

# Viewer.add_items(*AssemblyData.get_parts())
all_parts += AssemblyData.get_parts()

# STRUTS ----------------------------------------------------------------------
crv = strut.extract_curve((0.5, 0.), (0.5, 1.))
Viewer.add_entity(crv)

crv = jury.extract_curve((0.5, 0.), (0.5, 1.))
Viewer.add_entity(crv)

crv = other_strut.extract_curve((0.5, 0.), (0.5, 1.))
Viewer.add_entity(crv)

crv = other_jury.extract_curve((0.5, 0.), (0.5, 1.))
Viewer.add_entity(crv)

# VIEW ------------------------------------------------------------------------
compound = ShapeTools.make_compound(all_parts)
faces = ShapeTools.get_faces(compound)
for f in faces:
    Viewer.add_entity(f, 'random')
Viewer.show(False)

# MESH ------------------------------------------------------------------------
shape_to_mesh = ShapeTools.make_compound(all_parts)
MeshData.create_mesh('tbw mesh', shape_to_mesh)
MeshData.hypotheses.create_netgen_simple_2d('netgen hypo', 4., True)
MeshData.hypotheses.create_netgen_algo_2d('netgen algo')
MeshData.add_hypothesis('netgen hypo')
MeshData.add_hypothesis('netgen algo')

# Compute the mesh
mesh_start = time.time()
print('Computing mesh...')
status = MeshData.compute_mesh()
if not status:
    print('Failed to compute mesh')
else:
    print('Meshing complete in ', time.time() - mesh_start, ' seconds.')

Viewer.add_meshes(MeshData.get_active())
Viewer.show()
