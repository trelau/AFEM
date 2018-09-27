import time

from afem.graphics import Viewer
from afem.smesh import MeshGen, NetgenAlgo2D3D, NetgenHypothesis, MeshEditor
from afem.topology import BoxBySize

# Create a simple solid box
box = BoxBySize(50, 10, 10).solid

# Initialize the mesh generator
gen = MeshGen()

# Create a new mesh using the box as the shape
mesh = gen.create_mesh(box)

# 3-D volume controls
alg3d = NetgenAlgo2D3D(gen)
hyp3d = NetgenHypothesis(gen, max_size=0.5, min_size=0.5)

# Apply mesh controls
mesh.add_hypotheses([alg3d, hyp3d], box)

# Compute the mesh
print('Computing the mesh...')
start = time.time()
gen.compute(mesh)
t = time.time() - start
print('done in {} seconds.'.format(t))
print('Mesh elements: {}'.format(mesh.num_volumes))
print('Mesh nodes: {}'.format(mesh.num_nodes))

# View the mesh
gui = Viewer()
gui.add(mesh)
gui.start()
gui.clear()

# Upgrade mesh to quadratic
editor = MeshEditor(mesh)
print('Converting to quadratic...')
start = time.time()
status = editor.convert_to_quadratic()
t = time.time() - start
print('done in {} seconds.'.format(t))
print('Mesh elements: {}'.format(mesh.num_volumes))
print('Mesh nodes: {}'.format(mesh.num_nodes))
