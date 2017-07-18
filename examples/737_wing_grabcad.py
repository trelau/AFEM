from afem.fem import MeshData
from afem.geometry import CreateGeom, ProjectGeom
from afem.graphics import Viewer
from afem.io import StepImport
from afem.oml.wing import Wing
from afem.structure import AssemblyData, CreatePart, PartTools
from afem.topology import ShapeTools

# Import wing solid from STEP file
step_read = StepImport()
step_read.read(r'../models/boeing_737_wing_grabcad.step')
shape = step_read.shape
solids = ShapeTools.get_solids(shape)
# There should only be one solid that is the wing
solid = solids[0]

# Manually create a wing reference surface. Points came from examination in
# CAD program.
p1 = [0, 0, 0.41811]
p2 = [-299.374016, 0, -4.790998]
c1 = CreateGeom.linear_curve([p1, p2])

p1 = [-147.834646, 216.535433, 19.148031]
p2 = [-297.834646, 216.535433, 17.768032]
c2 = CreateGeom.linear_curve([p1, p2])

p1 = [-370.472441, 669.291339, 58.555188]
p2 = [-424.050551, 669.291339, 58.597981]
c3 = CreateGeom.linear_curve([p1, p2])

sref = CreateGeom.linear_surface([c1, c2, c3])
wing = Wing(solid)
wing.set_sref(sref)

# Root rib and fuselage radius
r = 148. / 2.
pln = CreateGeom.plane_by_axes([0, r, 0], 'xz')
root_rib = CreatePart.rib.by_sref('root rib', wing, pln)

# Tip rib
p0 = [-370.472441, 669.291339, 58.555188]
pln = CreateGeom.plane_by_axes(p0, 'xz')
tip_rib = CreatePart.rib.by_sref('tip rib', wing, pln)

# Rear spar
p1 = root_rib.eval_dx(0.7, is_local=True)
p2 = tip_rib.eval_dx(0.7, is_local=True)
rspar = CreatePart.spar.by_points('rear spar', wing, p1, p2)

# Front inboard spar
u = root_rib.local_to_global_u(0.15)
p1 = root_rib.eval_cref(u)
v = wing.vknots[1]
p2 = wing.eval(0.15, v)
inbd_fspar = CreatePart.spar.by_points('inbd front spar', wing, p1, p2)

# Front outboard spar
p1 = inbd_fspar.p2
u = tip_rib.local_to_global_u(0.15)
p2 = tip_rib.eval_cref(u)
outbd_fspar = CreatePart.spar.by_points('outbd front spar', wing, p1, p2)

# Kink rib
p1 = inbd_fspar.p2
p2 = inbd_fspar.p2
rspar.points_to_cref([p2])
pln = ShapeTools.plane_from_section(inbd_fspar, outbd_fspar, p2)
kink_rib = CreatePart.rib.by_points('kink rib', wing, p1, p2, pln)

# Inboard ribs
u2 = rspar.invert_cref(kink_rib.p2)
CreatePart.rib.along_curve('inbd rib', wing, rspar.cref,
                           inbd_fspar.sref, rspar.sref, 24., s1=12,
                           s2=-24, u2=u2)

# Outboard ribs
CreatePart.rib.along_curve('outbd rib', wing, rspar.cref,
                           outbd_fspar.sref, rspar.sref, 24., s1=24,
                           s2=-24, u1=u2)

# Front center spar
xz_pln = CreateGeom.plane_by_axes(axes='xz')
p2 = inbd_fspar.p1
p1 = inbd_fspar.p1
ProjectGeom.point_to_geom(p1, xz_pln, True)
pln = ShapeTools.plane_from_section(inbd_fspar, root_rib, p1)
fcspar = CreatePart.spar.by_points('front center spar', wing, p1, p2, pln)

# Rear center spar
p2 = rspar.p1
p1 = rspar.p1
ProjectGeom.point_to_geom(p1, xz_pln, True)
pln = ShapeTools.plane_from_section(rspar, root_rib, p1)
rcspar = CreatePart.spar.by_points('rear center spar', wing, p1, p2, pln)

# Center ribs
pac = CreateGeom.points_along_curve(rcspar.cref, npts=4)
i = 1
for p in pac.pnts[1:-1]:
    pln = CreateGeom.plane_by_axes(p, 'xz')
    name = ' '.join(['center rib', str(i)])
    CreatePart.rib.between_geom(name, wing, fcspar.sref,
                                rcspar.sref, pln)
    i += 1

wing_parts = AssemblyData.get_parts()

# Fuse wing parts and discard faces
PartTools.fuse_wing_parts(wing_parts)
PartTools.discard_faces(wing_parts)

# Skin
skin = CreatePart.skin.from_body('skin', wing)
skin.fuse(*wing_parts)

hs = ShapeTools.make_halfspace(tip_rib.sref, [0, 1e6, 0])

skin.discard(hs)

# View geometry
skin.set_transparency(0.5)
Viewer.add(*AssemblyData.get_parts())
Viewer.show(False)

# Mesh
print('Meshing the shape...')
the_shape = AssemblyData.prepare_shape_to_mesh()
MeshData.create_mesh('the mesh', the_shape)
MeshData.hypotheses.create_netgen_simple_2d('netgen', 4.)
MeshData.hypotheses.create_netgen_algo_2d('netgen algo')
MeshData.add_hypothesis('netgen')
MeshData.add_hypothesis('netgen algo')
MeshData.compute_mesh()

Viewer.add_meshes(MeshData.get_active())
Viewer.show()
