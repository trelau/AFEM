import time

from afem.exchange import brep
from afem.graphics import Viewer
from afem.smesh import *

fn = '../models/wing_body.brep'
shape = brep.read_brep(fn)

# MeshGems
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(shape)
alg2d = MeshGemsAlgo2D(the_gen)
hyp2d = MeshGemsHypo2D(the_gen, 4.)
the_mesh.add_hypotheses([alg2d, hyp2d])
print('Computing mesh with MeshGems...')
start = time.time()
the_gen.compute(the_mesh)
print('MeshGems complete in ', time.time() - start, ' seconds.')

gui = Viewer()
gui.add(the_mesh)
gui.start()
