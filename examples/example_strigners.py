from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart, PartTools
from asap.topology import ShapeTools

# Inputs
height = 3.
pitch = 8.

# Import model
fn = r'..\models\777-200LR.stp'
ImportVSP.step_file(fn)

# Get components.
wing = ImportVSP.get_body('Wing')
wing.set_transparency(0.5)

# Create some internal structure for skin boundary.
s1 = CreatePart.spar.by_parameters('spar1', wing, 0.15, 0.05, 0.15, 0.8)
s2 = CreatePart.spar.by_parameters('spar2', wing, 0.65, 0.05, 0.65, 0.8)
r1 = CreatePart.rib.by_points('rib1', wing, s1.p1, s2.p1)
r2 = CreatePart.rib.by_points('rib2', wing, s1.p2, s2.p2)

PartTools.fuse_wing_parts([s1, s2, r1, r2])
PartTools.discard_faces([s1, s2, r1, r2])

skin = CreatePart.skin.from_body('skin', wing)
skin.fuse(s1, s2, r1, r2)
skin.discard(wing.sref)
skin.fix()

# Get one of the skin shells and find a common edge to distribute points
# along.
upper_shell = ShapeTools.get_shells(skin)[1]
skin = CreatePart.surface_part('skin', upper_shell)
edge = skin.shared_edges(r1)[0]
pnts = ShapeTools.points_along_edge(edge, pitch)[1:-1]

# Create stringers parallel to the rear spar.
vn = s2.sref.norm(0, 0)
for p in pnts:
    pln = CreateGeom.plane_by_normal(p, vn)
    stringer = CreatePart.stringer.by_section('stringer', skin, pln, height)
    # Some stringers needs fixed to resolve topology.
    stringer.fix()
    Viewer.add_items(stringer)

Viewer.add_items(skin)
Viewer.show()
