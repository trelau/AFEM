import time

from afem.config import Settings
from afem.fem import MeshAPI
from afem.geometry import *
from afem.graphics.display import display_shape
from afem.io import ImportVSP
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

fname = r'..\models\777-200LR.stp'
ImportVSP.step_file(fname)
wing = ImportVSP.get_body('Wing')
fuselage = ImportVSP.get_body('Fuselage')
htail = ImportVSP.get_body('Htail')

wing.set_transparency(0.5)
fuselage.set_transparency(0.5)

# WING
wing_assy = AssemblyAPI.create_assy('wing assy')

rib_assy = AssemblyAPI.create_assy('rib assy', wing_assy)

root = RibByParameters('root', 0.15, 0.05, 0.70, 0.05, wing).rib
tip = RibByParameters('tip', 0.15, 0.2, 0.70, 0.2, wing).rib

spar_assy = AssemblyAPI.create_assy('spar assy', wing_assy)

fspar = SparByPoints('fspar', root.p1, tip.p1, wing).spar
rspar = SparByPoints('rspar', root.p2, tip.p2, wing).spar

rib_assy.activate()
RibsAlongCurveByNumber('rib', rspar.cref, 5, fspar, rspar, wing, d1=30.,
                       d2=-50)

# Use FuseAssemblies
FuseAssemblies([rib_assy, spar_assy])
DiscardByCref(wing_assy.get_parts())

wing_assy.activate()
skin = SkinByBody('wing skin', wing).skin
parts = rib_assy.get_parts() + spar_assy.get_parts()
skin.fuse(*parts)
skin.discard_by_dmin(wing.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

# FUSELAGE
fuse_assy = AssemblyAPI.create_assy('fuse assy')

skin = SkinByBody('fuselage skin', fuselage).skin
skin.set_transparency(0.5)

pln = PlaneByAxes(root.point_from_parameter(0.5, is_rel=True), 'yz').plane
frame = FrameByPlane('frame', pln, fuselage, 3.).frame

# Causing issues with fusing assy
# vol = VolumesFromShapes(wing_assy.get_parts()).shape
# skin.cut(vol)
# frame.cut(vol)

skin.fuse(frame)

# HTAIL
htail_assy = AssemblyAPI.create_assy('htail assy')

root = RibByParameters('root', 0.15, 0.05, 0.70, 0.05, htail).rib
tip = RibByParameters('tip', 0.15, 0.5, 0.70, 0.5, htail).rib

fspar = SparByPoints('fspar', root.p1, tip.p1, htail).spar
rspar = SparByPoints('rspar', root.p2, tip.p2, htail).spar

pln = PlaneByAxes().plane
RibsAlongCurveByNumber('rib', rspar.cref, 5, fspar, rspar, htail, d1=12,
                       d2=-12, ref_pln=pln)

parts = htail_assy.get_parts()
FuseSurfacePartsByCref(parts)
DiscardByCref(parts)

skin = SkinByBody('htail skin', htail).skin
skin.fuse(*parts)
skin.discard_by_dmin(htail.sref_shape, 1.)
skin.fix()
skin.set_transparency(0.5)

# JOIN
start = time.time()
bop = FuseAssemblies([wing_assy, fuse_assy, htail_assy])
print(time.time() - start)

shape = bop.fused_shape
tool = ExploreFreeEdges(shape)
print(tool.free_edges)

parts = AssemblyAPI.get_master().get_parts()

# MESH
master_assy = AssemblyAPI.get_master()
shape_to_mesh = master_assy.prepare_shape_to_mesh()
the_mesh = MeshAPI.create_mesh('wing-box mesh', shape_to_mesh)

# Use a single global hypothesis based on local length.
MeshAPI.hypotheses.create_netgen_simple_2d('netgen hypo', 4.)
MeshAPI.hypotheses.create_netgen_algo_2d('netgen algo')
MeshAPI.add_hypothesis('netgen hypo')
MeshAPI.add_hypothesis('netgen algo')

# Compute the mesh
mesh_start = time.time()
print('Computing mesh...')
status = MeshAPI.compute_mesh()
if not status:
    print('Failed to compute mesh')
else:
    print('Meshing complete in ', time.time() - mesh_start, ' seconds.')

display_shape(None, the_mesh.object, *tool.free_edges)
