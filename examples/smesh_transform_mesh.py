from OCCT.gp import gp_Trsf

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

# Mesh using composite side algorithm to avoid making a vertex at edge
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(face)
hyp1d = LocalLength1D(the_gen, 1)
alg1d = Regular1D(the_gen)
the_mesh.add_hypotheses([hyp1d, alg1d], wire)
hyp2d = QuadrangleHypo2D(the_gen)
alg2d = QuadrangleAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp2d, alg2d], face)

the_gen.compute(the_mesh, face)

editor = MeshEditor(the_mesh)

# TODO Transate a mesh
pz = Point(0, 0, 10)
trsf = gp_Trsf()
trsf.SetTranslation(p1, pz)

new_mesh = the_gen.create_mesh(face)

editor.transform(trsf, copy=True)

gui = Viewer()
gui.add(the_mesh)
gui.start()
