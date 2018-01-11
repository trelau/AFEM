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
wire = WireByPoints([p1, p2, p3, p4, p5], True).wire
face = FaceByPlanarWire(wire).face

# Mesh using composite side algorithm to avoid making a vertex at edge
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(face)
hyp1d = LocalLength1D(the_gen, 4)
alg1d = CompositeSide1D(the_gen)
the_mesh.add_hypotheses([hyp1d, alg1d], wire)
hyp2d = QuadrangleHypo2D(the_gen)
alg2d = QuadrangleAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp2d, alg2d], face)

# Add enforced node in middle
hyp2d.set_enforced_nodes([face], [(5, 5, 0)])

the_gen.compute(the_mesh, face)

v = Viewer()
v.view.view_top()
for vert in ExploreShape.get_vertices(face):
    v.add(vert)
v.view.display_mesh(the_mesh.object, 2)
v.start()
