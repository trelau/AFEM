from OCC.GC import GC_MakeCircle

from afem.geometry import CreateGeom
from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.structure import CreatePart, PartTools
from afem.topology import ShapeTools

# Inputs (inches)
height = 6.
r0 = 36.
fn = r'..\models\777-200LR.stp'

# Import model
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

# Get one of the skin shells.
upper_shell = ShapeTools.get_shells(skin)[1]
skin = CreatePart.surface_part('skin', upper_shell)

# Define circles of increasing radius at root of root rib to generate curved
#  stringers.
line = CreateGeom.line_by_points(r1.p1, r1.p2)
p0 = line.eval(4. * r1.cref.length)
vn = CreateGeom.direction_by_xyz(0, 0, 1)

v = Viewer()
# Takes a while to run...
for i in range(22, 31):
    geom_circle = GC_MakeCircle(p0, vn, r0 * (i + 1)).Value().GetObject()
    basis_shape = ShapeTools.make_prism(geom_circle, [0, 0, 200])
    stringers = CreatePart.stringer.by_sections('stringer', skin, basis_shape,
                                                height, 30.)
    for stringer in stringers:
        print(stringer, ShapeTools.is_valid(stringer))
        if not ShapeTools.is_valid(stringer):
            stringer.fix()
        v.add(stringer)

v.add(skin)
v.set_display_shapes()
v.show()
