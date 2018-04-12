from afem.geometry import *
from afem.graphics import *
from afem.topology import *

gui = Viewer()

# Create a box by size
builder = BoxBySize(10, 10, 10)
box = builder.solid
box.set_transparency(0.5)

# Create a cylinder partially inside the box
circle = CircleByNormal((5, 5, 5), (0, 0, 1), 2).circle
face = FaceByPlanarWire(circle).face
cyl = SolidByDrag(face, (0, 0, 15)).solid

# View the two shapes
gui.add(box, cyl)
gui.start()
gui.clear()

# Fuse the shapes
fuse = FuseShapes(box, cyl)
fused_shape = fuse.shape
fused_shape.set_transparency(0.5)

gui.add(fused_shape)
gui.start()
gui.clear()

# Cut the cylinder from the box
cut = CutShapes(box, cyl)
cut_shape = cut.shape

gui.add(cut_shape)
gui.start()
gui.clear()

# Common material between the two shapes
common = CommonShapes(box, cyl)
common_shape = common.shape

# Show original box for reference
gui.add(common_shape, box)
gui.start()
gui.clear()

# Intersect the shapes
sec = IntersectShapes(box, cyl)
sec_shape = sec.shape

# Original shapes shown for reference
gui.add(sec_shape, box, cyl)
gui.start()
gui.clear()

# Split the box with the cylinder. The resulting shape is a compound with two
# solids.
split = SplitShapes(box, cyl)
split_shape = split.shape
split_shape.set_transparency(0.5)

gui.add(split_shape)
gui.start()
gui.clear()

# Locally split one face of the box with a plane
pln = PlaneByAxes((5, 5, 5), 'xz').plane

local = LocalSplit(builder.front_face, pln, box)
local_shape = local.shape
local_shape.set_transparency(0.5)

gui.add(local_shape)
gui.start()
gui.clear()

# Offset the box
offset = OffsetShape(box, 2)
offset_shape = offset.shape
offset_shape.set_transparency(0.5)

gui.add(box, offset_shape)
gui.start()
gui.clear()

# Rebuild the box with the results of the cut tool.
rebuild = RebuildShapeByTool(box, cut)
new_shape = rebuild.new_shape

gui.add(new_shape)
gui.start()

# Check the new shape for errors
check = CheckShape(new_shape)
print('Shape is valid:', check.is_valid)
print('Shape type:', new_shape.shape_type)

# Since a face is removed it is no longer a valid solid but a shell. Try to
# fix the shape.
fix = FixShape(new_shape)
fixed_shape = fix.shape

check = CheckShape(fixed_shape)
print('Shape is valid:', check.is_valid)
print('Shape type:', fixed_shape.shape_type)

gui.add(fixed_shape)
gui.start()

# Find free edges of a shape
tool = ExploreFreeEdges(fixed_shape)

gui.add(*tool.free_edges)
gui.start()
