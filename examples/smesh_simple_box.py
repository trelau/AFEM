from OCCT.BRepPrimAPI import BRepPrimAPI_MakeBox

from afem.graphics import Viewer
from afem.smesh import (MeshGen, NetgenAlgo2D, NetgenSimple2D, MaxLength1D,
                        Regular1D)
from afem.topology import ExploreShape

# Create a simple solid box
box = BRepPrimAPI_MakeBox(10, 10, 10).Solid()

# Get a list of faces and edges of the shape for later use
faces = ExploreShape.get_faces(box)
edges = ExploreShape.get_edges(box)

# Initialize the mesh generator
gen = MeshGen()

# Create a new mesh using the box as the shape
mesh = gen.create_mesh(box)

# Define algorithms and hypotheses
alg2d = NetgenAlgo2D(gen)
hyp2d_1 = NetgenSimple2D(gen, 1.)
hyp2d_2 = NetgenSimple2D(gen, 1., allow_quads=False)
alg1d = Regular1D(gen)
hyp1d = MaxLength1D(gen, 0.25)

# Add them to the mesh
mesh.add_hypotheses([alg2d, hyp2d_1])
mesh.add_hypothesis(hyp2d_2, faces[-1])
mesh.add_hypotheses([alg1d, hyp1d], edges[-1])

# Compute the mesh
gen.compute(mesh)

# View the mesh
v = Viewer()
v.add(mesh)
v.start()
v.clear()

# Get sub-meshes from sub-shapes
face_submesh = mesh.get_submesh(faces[0])
edge_submesh = mesh.get_submesh(edges[0])

# View the face sub-mesh (2-D elements)
v.add(face_submesh)
v.start()
v.clear()

# View the edge sub-mesh (1-D elements)
v.add(edge_submesh)
v.start()
