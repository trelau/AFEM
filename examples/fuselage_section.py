from OCC.BRepPrimAPI import BRepPrimAPI_MakeCylinder

from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.oml.fuselage import Fuselage
from asap.structure import CreatePart, PartTools, AssemblyData
from asap.topology import ShapeTools
from asap.fem import MeshData

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
cylinder = BRepPrimAPI_MakeCylinder(radius, length).Shape()
fuselage = Fuselage(cylinder)

# Skin
skin = CreatePart.skin.from_body('skin', fuselage)

# Trim off closed ends of skin since it came from a solid cylinder.
pln1 = CreateGeom.plane_by_axes([0, 0, 0], 'xy')
box = ShapeTools.box_from_plane(pln1, 1e6, 1e6, -1e6)
skin.cut(box)

pln2 = CreateGeom.plane_by_axes([0, 0, length], 'xy')
box = ShapeTools.box_from_plane(pln2, 1e6, 1e6, 1e6)
skin.cut(box)

# Floor
pln = CreateGeom.plane_by_axes([0, main_floor_yloc, 0], 'xz')
main_floor = CreatePart.floor.by_sref('main floor', fuselage, pln)

pln = CreateGeom.plane_by_axes([0, cargo_floor_yloc, 0], 'xz')
cargo_floor = CreatePart.floor.by_sref('cargo floor', fuselage, pln)

# Frames
frames = CreatePart.frame.between_planes('frame', fuselage, [pln1, pln2],
                                         frame_height, frame_spacing)

# Floor beams and posts
rev_cylinder = cylinder.Reversed()
above_floor = ShapeTools.make_prism(main_floor, [0, 2 * diameter, 0])
below_cargo_floor = ShapeTools.make_prism(cargo_floor, [0, -60, 0])
pln1 = CreateGeom.plane_by_axes([-.667 * radius, 0, 0], 'yz')
face1 = ShapeTools.face_from_plane(pln1, -diameter, diameter, 0, length)
pln2 = CreateGeom.plane_by_axes([.667 * radius, 0, 0], 'yz')
face2 = ShapeTools.face_from_plane(pln2, -diameter, diameter, 0, length)
i = 1
for frame in frames:
    # Beam
    shape = ShapeTools.bsection(main_floor, frame.sref)
    shape = ShapeTools.make_prism(shape, [0, -floor_beam_height, 0])
    name = ' '.join(['floor beam', str(i)])
    beam = CreatePart.surface_part(name, shape)
    beam.cut(rev_cylinder)
    # Post
    name = ' '.join(['left floor post', str(i)])
    post = CreatePart.curve_part(name, shape1=face1, shape2=frame.sref)
    post.cut(above_floor)
    post.cut(rev_cylinder)
    name = ' '.join(['right floor post', str(i)])
    post = CreatePart.curve_part(name, shape1=face2, shape2=frame.sref)
    post.cut(above_floor)
    post.cut(rev_cylinder)
    # Create segment beneath cargo floor and merge with frame.
    shape = ShapeTools.bcommon(below_cargo_floor, frame.sref)
    shape = ShapeTools.bcut(shape, rev_cylinder)
    frame.merge(shape, True)
    i += 1

main_floor.set_transparency(0.5)
cargo_floor.set_transparency(0.5)

# Cut the main floor with post planes.
main_floor.cut(pln1)
main_floor.cut(pln2)

# Fuse all parts together.
PartTools.split_parts(AssemblyData.get_parts())

Viewer.add_items(*AssemblyData.get_parts())
Viewer.show(False)

# Mesh
the_shape = AssemblyData.prepare_shape_to_mesh()
MeshData.create_mesh('the mesh', the_shape)
MeshData.hypotheses.create_netgen_simple_2d('netgen', 4.)
MeshData.hypotheses.create_netgen_algo_2d('netgen algo')
MeshData.add_hypothesis('netgen')
MeshData.add_hypothesis('netgen algo')
# edges = ShapeTools.get_edges(the_shape, True)
# MeshData.hypotheses.create_local_length_1d('local length', 4)
# MeshData.hypotheses.create_regular_1d('algo 1d')
# MeshData.add_hypothesis('local length', edges)
# MeshData.add_hypothesis('algo 1d', edges)
MeshData.compute_mesh()

Viewer.add_meshes(MeshData.get_active())
Viewer.show()
