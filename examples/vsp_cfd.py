"""
Demonstrate generating a watertight solid from OpenVSP models for CFD purposes.
"""
import time

from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer
from afem.smesh import *
from afem.topology import FixShape, FuseShapes

Settings.log_to_console()

# Import an OpenVSP STEP file. If generated using the modified version that
# includes metadata, each Body will be retrievable by its component name.
fn = '../models/777-200LR.stp'

# Import the OpenVSP STEP file
vsp_import = ImportVSP(fn)

# View the bodies
gui = Viewer()
solids = []
for body in vsp_import.all_bodies:
    gui.add(body)
    solids.append(body.shape)
gui.start()
gui.clear()

# The Boolean fuse operation must have arguments and tools. Pick one solid to
# be an argument and the others tools. So far the choice of argument and tools
# has not impacted the final fused shape. If you're able to pick the arguments
# and tools, perhaps picking the largest/most connected solid is best (like
# a fuselage)?
args = [solids[0]]
tools = solids[1:]

# Fuse the shapes
fuse = FuseShapes()
fuse.set_args(args)
fuse.set_tools(tools)

print('Fusing the shapes...')
start = time.time()
fuse.build()
print('Complete in ', time.time() - start, 'seconds.')

# Get the resulting shape
shape = fuse.shape

# Inspect the average shape tolerance and limit for robustness. This is not
# a necessary step, especially if modeling operations are complete.
print('Tol after fuse: ', shape.tol_avg)
FixShape.limit_tolerance(shape)
print('Tol after limiting: ', shape.tol_avg)
shape.set_transparency(0.5)

# Display the fused shape
gui.add(shape)
gui.start()
gui.clear()

# Mesh the shape with Netgen
gen = MeshGen()
mesh = gen.create_mesh(shape)

alg = NetgenAlgo2D(gen)
hyp = NetgenHypo2D(gen, max_size=10., min_size=1., surface_curvature=True)

mesh.add_hypotheses([alg, hyp])

print('Meshing the shape...')
start = time.time()
gen.compute(mesh)
print('Complete in ', time.time() - start, 'seconds.')

# Native export options from SMESH include CGNS, DAT, GMSH, STL, UNV, SAUV
# See "mesh.object.Export*" methods which can be exposed to AFEM by request.

# Shape export include BREP, STEP, IGES, and STL and can be found in the
# afem.exchange sub-package.

# View the shape using flat shading. Need to figure out how to better render
# meshes.
gui.add(mesh)
gui.start()
