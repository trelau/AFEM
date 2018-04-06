from afem.geometry import *
from afem.graphics import Viewer
from afem.smesh import *
from afem.topology import *

# Create a face with a vertex in middle of one edge
p1 = Point()
p2 = Point(10, 0, 0)
p3 = Point(10, 10, 0)
p4 = Point(5, 10, 0)
p5 = Point(0, 10, 0)
wire = Wire.by_points([p1, p2, p3, p4, p5], True)
face = Face.by_wire(wire)

# Mesh using composite side algorithm to avoid making a vertex at edge
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(face)
hyp1d = LocalLength1D(the_gen, 4)
alg1d = CompositeSide1D(the_gen)
the_mesh.add_hypotheses([hyp1d, alg1d], wire)
hyp2d = NetgenSimple2D(the_gen, 1)
alg2d = NetgenAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp2d, alg2d], face)
the_gen.compute(the_mesh, face)

# Get the composite side
for e in wire.edges:
    fside = alg1d.get_face_side(the_mesh, e, face)
    if fside.num_nodes:
        break

gui = Viewer()
gui.view_top()
for vert in face.vertices:
    gui.add(vert)
gui.add(the_mesh)
gui.start()
