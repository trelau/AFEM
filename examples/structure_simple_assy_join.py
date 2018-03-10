import time

from afem.config import Settings
from afem.exchange import ImportVSP
from afem.geometry import *
from afem.graphics import Viewer
from afem.smesh import *
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

fname = r'..\models\777-200LR.stp'
vsp_import = ImportVSP(fname)
wing = vsp_import.get_body('Wing')
fuselage = vsp_import.get_body('Fuselage')
htail = vsp_import.get_body('Htail')

wing.set_transparency(0.5)
fuselage.set_transparency(0.5)

# WING
wing_group = GroupAPI.create_group('wing group')

rib_group = GroupAPI.create_group('rib group', wing_group)

root = RibByParameters('root', 0.15, 0.05, 0.70, 0.05, wing).rib
tip = RibByParameters('tip', 0.15, 0.2, 0.70, 0.2, wing).rib

spar_group = GroupAPI.create_group('spar group', wing_group)

fspar = SparByPoints('fspar', root.p1, tip.p1, wing).spar
rspar = SparByPoints('rspar', root.p2, tip.p2, wing).spar

rib_group.activate()
RibsAlongCurveByNumber('rib', rspar.cref, 5, fspar.shape, rspar.shape, wing,
                       d1=30., d2=-50)

# Use FuseGroups
FuseGroups([rib_group, spar_group])
DiscardByCref(wing_group.get_parts())

wing_group.activate()
skin = SkinByBody('wing skin', wing).skin
parts = rib_group.get_parts() + spar_group.get_parts()
skin.fuse(*parts)
skin.discard_by_dmin(wing.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

# FUSELAGE
fuse_group = GroupAPI.create_group('fuse group')

skin = SkinByBody('fuselage skin', fuselage).skin
skin.set_transparency(0.5)

pln = PlaneByAxes(root.point_from_parameter(0.5, is_rel=True), 'yz').plane
frame = FrameByPlane('frame', pln, fuselage, 3.).frame

# Causing issues with fusing group
# vol = VolumesFromShapes(wing_group.get_parts()).shape
# skin.cut(vol)
# frame.cut(vol)

skin.fuse(frame)

# HTAIL
htail_group = GroupAPI.create_group('htail group')

root = RibByParameters('root', 0.15, 0.05, 0.70, 0.05, htail).rib
tip = RibByParameters('tip', 0.15, 0.5, 0.70, 0.5, htail).rib

fspar = SparByPoints('fspar', root.p1, tip.p1, htail).spar
rspar = SparByPoints('rspar', root.p2, tip.p2, htail).spar

pln = PlaneByAxes().plane
RibsAlongCurveByNumber('rib', rspar.cref, 5, fspar.shape, rspar.shape, htail,
                       d1=12, d2=-12, ref_pln=pln)

parts = htail_group.get_parts()
FuseSurfacePartsByCref(parts)
DiscardByCref(parts)

skin = SkinByBody('htail skin', htail).skin
skin.fuse(*parts)
skin.discard_by_dmin(htail.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

# JOIN
start = time.time()
bop = FuseGroups([wing_group, fuse_group, htail_group])
print(time.time() - start)

shape = bop.shape
tool = ExploreFreeEdges(shape)
print(tool.free_edges)

parts = GroupAPI.get_master().get_parts()

# MESH
master_group = GroupAPI.get_master()
shape_to_mesh = master_group.prepare_shape_to_mesh()
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(shape_to_mesh)

# Use a single global hypothesis based on local length.
hyp = NetgenSimple2D(the_gen, 4.)
alg = NetgenAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp, alg], shape_to_mesh)

# Compute the mesh
mesh_start = time.time()
print('Computing mesh...')
status = the_gen.compute(the_mesh, shape_to_mesh)
if not status:
    print('Failed to compute mesh')
else:
    print('Meshing complete in ', time.time() - mesh_start, ' seconds.')


v = Viewer()
v.add(*tool.free_edges)
v.display_mesh(the_mesh.object, 2)
v.start()
