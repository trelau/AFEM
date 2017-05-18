import time

from asap.fem import MeshData
from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import AssemblyData, CreatePart, PartTools

# Inputs
fname = r'..\models\N2A_nosplit.stp'
fd_length = 10. * 12.
cabin_length = 100. * 12.

# Import model
ImportVSP.step_file(fname)
wing = ImportVSP.get_body('Wing_Body')
other_wing = ImportVSP.get_body('Wing_Body.1')
vtail = ImportVSP.get_body('Vertical_Tails')
other_vtail = ImportVSP.get_body('Vertical_Tails.2')
for name in ImportVSP.get_bodies():
    body = ImportVSP.get_body(name)
    body.set_transparency(0.5)
    body.set_color(0.5, 0.5, 0.5)
    Viewer.add_items(body)

# Centerbody structure
sref = CreateGeom.plane_by_axes([fd_length, 0, 0], 'yz')
fd_bh = CreatePart.bulkhead.by_sref('flight deck bulkhead', wing, sref)

sref = CreateGeom.plane_by_axes([fd_length + cabin_length, 0, 0], 'yz')
rear_cabin_bh = CreatePart.bulkhead.by_sref('rear cabin bulkhead', wing, sref)

# Outboard wing structure

# Root rib
outbd_root_rib = CreatePart.rib.by_parameters('outbd root rib', wing,
                                              0.15, 0.5, 0.65, 0.5)

# Tip rib
tip_rib = CreatePart.rib.by_parameters('tip rib', wing,
                                       0.15, 0.995, 0.65, 0.995)

# Outbd front spar
outbd_fspar = CreatePart.spar.by_points('outbd front spar', wing,
                                        outbd_root_rib.p1, tip_rib.p1)

# Outbd rear spar
outbd_rspar = CreatePart.spar.by_points('outbd rear spar', wing,
                                        outbd_root_rib.p2, tip_rib.p2)

# Outbd ribs
outbd_ribs = CreatePart.rib.along_curve('outbd rib', wing, outbd_rspar.cref,
                                        outbd_fspar.sref, outbd_rspar.sref,
                                        72., s1=30., s2=-30.)

internal_parts = AssemblyData.get_parts()

PartTools.fuse_wing_parts(internal_parts)
PartTools.discard_faces(internal_parts)

Viewer.add_items(*internal_parts)

Viewer.show(False)

# MESH --------------------------------------------------------------------
# Initialize
shape_to_mesh = AssemblyData.prepare_shape_to_mesh()
MeshData.create_mesh('N2A mesh', shape_to_mesh)

# Use a single global hypothesis based on local length.
MeshData.hypotheses.create_netgen_simple_2d('netgen hypo', 4.)
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
