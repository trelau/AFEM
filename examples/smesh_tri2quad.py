from afem.graphics import Viewer
from afem.smesh import *
from afem.topology import *

box = BoxBySize(10, 10, 10).solid

the_gen = MeshGen()
the_mesh = the_gen.create_mesh(box)
hyp2d = NetgenSimple2D(the_gen, 1, allow_quads=False)
alg2d = NetgenAlgo2D(the_gen)
the_mesh.add_hypotheses([hyp2d, alg2d])
the_gen.compute(the_mesh)

editor = MeshEditor(the_mesh)
editor.tri_to_quad()

editor.smooth(iters=20)

gui = Viewer()
gui.add(the_mesh)
gui.start()
