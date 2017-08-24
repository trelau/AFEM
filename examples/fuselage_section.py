from afem.fem import MeshAPI
from afem.geometry import *
from afem.graphics import Viewer
from afem.oml import *
from afem.structure import *
from afem.topology import *

# Inputs
diameter = 244
length = 360
main_floor_yloc = -12
cargo_floor_yloc = -108
frame_height = 3.5
frame_spacing = 24
floor_beam_height = 6

# Calculate
radius = diameter / 2.

# Create a solid cylinder to represent fuselage section.
cylinder = SolidByCylinder(radius, length).solid
fuselage = Body(cylinder)

# Skin
skin = SkinByBody('skin', fuselage).skin

# Trim off closed ends of skin since it came from a solid cylinder.
pln1 = PlaneByAxes(axes='xy').plane
box = SolidByPlane(pln1, 1e6, 1e6, -1e6).solid
skin.cut(box)

pln2 = PlaneByAxes((0., 0., length), 'xy').plane
box = SolidByPlane(pln2, 1e6, 1e6, 1e6).solid
skin.cut(box)

# Floor
pln = PlaneByAxes((0., main_floor_yloc, 0.), 'xz').plane
main_floor = FloorBySurface('main floor', pln, fuselage).floor

pln = PlaneByAxes((0., cargo_floor_yloc, 0.), 'xz').plane
cargo_floor = FloorBySurface('cargo floor', pln, fuselage).floor

# Frames
frames = FramesBetweenPlanesByDistance('frame', pln1, pln2, frame_spacing,
                                       fuselage, frame_height).frames

# Floor beams and posts
rev_cylinder = cylinder.Reversed()
above_floor = ShapeByDrag(main_floor, (0., 2. * diameter, 0.)).shape
below_cargo_floor = ShapeByDrag(cargo_floor, (0., -60., 0.)).shape

pln1 = PlaneByAxes((-.667 * radius, 0., 0.), 'yz').plane
face1 = FaceByPlane(pln1, 0., length, -diameter, diameter).face

pln2 = PlaneByAxes((.667 * radius, 0., 0.), 'yz').plane
face2 = FaceByPlane(pln2, 0., length, -diameter, diameter).face

i = 1
for frame in frames:
    # Beam
    shape = IntersectShapes(main_floor, frame.sref).shape
    shape = ShapeByDrag(shape, (0., -floor_beam_height, 0.)).shape
    name = ' '.join(['floor beam', str(i)])
    beam = SurfacePart(name, shape)
    beam.cut(rev_cylinder)

    # Post
    name = ' '.join(['left floor post', str(i)])
    shape = IntersectShapes(face1, frame.sref).shape
    post = CurvePart(name, shape)
    post.cut(above_floor)
    post.cut(rev_cylinder)

    name = ' '.join(['right floor post', str(i)])
    shape = IntersectShapes(face2, frame.sref).shape
    post = CurvePart(name, shape)
    post.cut(above_floor)
    post.cut(rev_cylinder)

    # Create segment beneath cargo floor and merge with frame.
    frame_pln_face = FaceBySurface(frame.sref).face
    shape = CommonShapes(below_cargo_floor, frame_pln_face).shape
    shape = CutShapes(shape, rev_cylinder).shape
    frame.merge(shape, True)
    i += 1

main_floor.set_transparency(0.5)
cargo_floor.set_transparency(0.5)

all_parts = AssemblyAPI.get_parts(order=True)

# Split all parts together
join = SplitParts(all_parts)

Viewer.add(*all_parts)
Viewer.show(False)

# Mesh
the_shape = AssemblyAPI.prepare_shape_to_mesh()
MeshAPI.create_mesh('the mesh', the_shape)
MeshAPI.hypotheses.create_netgen_simple_2d('netgen', 4.)
MeshAPI.hypotheses.create_netgen_algo_2d('netgen algo')
MeshAPI.add_hypothesis('netgen')
MeshAPI.add_hypothesis('netgen algo')
# edges = ShapeTools.get_edges(the_shape, True)
# MeshData.hypotheses.create_local_length_1d('local length', 4)
# MeshData.hypotheses.create_regular_1d('algo 1d')
# MeshData.add_hypothesis('local length', edges)
# MeshData.add_hypothesis('algo 1d', edges)
MeshAPI.compute_mesh()

Viewer.add_meshes(MeshAPI.get_active())
Viewer.show()
