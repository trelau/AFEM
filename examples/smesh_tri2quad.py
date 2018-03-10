from OCCT.BRepPrimAPI import BRepPrimAPI_MakeBox

from afem.graphics import Viewer
from afem.smesh import *

box = BRepPrimAPI_MakeBox(10, 10, 10).Solid()

the_gen = MeshGen()
the_mesh = the_gen.create_mesh(box)
hyp2d = NetgenSimple2D(the_gen, 1, allow_quads=False)
alg2d = NetgenAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp2d, alg2d])
the_gen.compute(the_mesh)

editor = MeshEditor(the_mesh)
editor.tri_to_quad()

editor.smooth(iters=20)

v = Viewer()
v.display_mesh(the_mesh.object, 2)
v.start()
