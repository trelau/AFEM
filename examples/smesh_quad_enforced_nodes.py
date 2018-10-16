from afem.geometry import *
from afem.graphics import Viewer
from afem.smesh import *
from afem.topology import *

# Create a face
p1 = Point()
p2 = Point(10, 0, 0)
p3 = Point(10, 10, 0)
p4 = Point(0, 10, 0)
wire = Wire.by_points([p1, p2, p3, p4], True)
face = Face.by_wire(wire)

# Mesh
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(face)
hyp1d = LocalLength1D(the_gen, 2)
alg1d = Regular1D(the_gen)
the_mesh.add_hypotheses([hyp1d, alg1d], wire)
hyp2d = QuadrangleHypo2D(the_gen)
alg2d = QuadrangleAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp2d, alg2d], face)

# Add enforced node in middle
p = Point(5, 5, 0)
hyp2d.set_enforced_nodes([face], [p])

the_gen.compute(the_mesh, face)

gui = Viewer()
gui.view_top()
gui.add(p)
gui.add(the_mesh)
gui.start()
