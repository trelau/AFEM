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
cylinder = CylinderByAxis(radius, length).solid
fuselage = Body(cylinder)

# Skin
skin = SkinByBody('skin', fuselage).part

# Trim off closed ends of skin since it came from a solid cylinder.
pln1 = PlaneByAxes(axes='xy').plane
box = SolidByPlane(pln1, 1e6, 1e6, -1e6).solid
skin.cut(box)

pln2 = PlaneByAxes((0., 0., length), 'xy').plane
box = SolidByPlane(pln2, 1e6, 1e6, 1e6).solid
skin.cut(box)

# Floor
pln = PlaneByAxes((0., main_floor_yloc, 0.), 'xz').plane
main_floor = FloorByShape('main floor', pln, fuselage).part

pln = PlaneByAxes((0., cargo_floor_yloc, 0.), 'xz').plane
cargo_floor = FloorByShape('cargo floor', pln, fuselage).part

# Frames
frames = FramesBetweenPlanesByDistance('frame', pln1, pln2, frame_spacing,
                                       fuselage, frame_height).parts

# Floor beams and posts
rev_cylinder = cylinder.reversed()
above_floor = ShapeByDrag(main_floor.shape, (0., 2. * diameter, 0.)).shape
below_cargo_floor = ShapeByDrag(cargo_floor.shape, (0., -60., 0.)).shape

pln1 = PlaneByAxes((-.667 * radius, 0., 0.), 'yz').plane
face1 = FaceByPlane(pln1, -diameter, diameter, 0., length).face

pln2 = PlaneByAxes((.667 * radius, 0., 0.), 'yz').plane
face2 = FaceByPlane(pln2, -diameter, diameter, 0., length).face

i = 1
for frame in frames:
    # Beam
    shape = IntersectShapes(main_floor.shape, frame.sref).shape
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

all_parts = GroupAPI.get_parts(order=True)

# Split all parts together
join = SplitParts(all_parts)

# Mesh
the_mesh = MeshVehicle(4.)

# Mapped quads applied to applicable faces
the_mesh.set_quadrangle_2d(the_mesh.shape)

print('Computing the mesh...')
the_mesh.compute()

# View
gui = Viewer()
gui.add(GroupAPI.get_master())
gui.start()
gui.clear()
gui.add(the_mesh)
gui.start()
