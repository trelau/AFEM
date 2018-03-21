from __future__ import print_function

from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.oml import Body
from afem.smesh import *
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Import model
fname = r'..\models\tbw.xbf'
bodies = Body.load_bodies(fname)

# Get the fuselage and define some bulkheads
fuselage = bodies['Fuselage']
fuselage.set_transparency(0.5)

pln1 = PlaneByAxes((100, 0, 0), 'yz').plane
pln2 = PlaneByAxes((400, 0, 0), 'yz').plane

planes = PlanesBetweenPlanesByNumber(pln1, pln2, 3).planes

bulkheads = []
for pln in planes:
    bh = BulkheadBySurface('bh', pln, fuselage).bulkhead
    bulkheads.append(bh)

# Define some edges to form a beam cross section. Build the section oriented
# with the started point and orientation in mind.
e1 = EdgeByPoints((100, -3, 4), (100, 3, 4)).edge
e2 = EdgeByPoints((100, -3, -4), (100, 3, -4)).edge
e3 = EdgeByPoints((100, 0, 4), (100, 0, -4)).edge
bop = FuseShapes()
bop.set_args([e1, e2])
bop.set_tools([e3])
bop.build()
profile = bop.shape

# Define the spine (i.e., path)
spine = EdgeByPoints((100, 0, 0), (400, 0, 0)).edge

# Build the beam
ibeam = Beam2DBySweep('beam', spine, profile).beam2d

# Fuse all the surface parts
FuseSurfaceParts(bulkheads, [ibeam])

v = Viewer()
v.add(fuselage, ibeam, *bulkheads)
print('Press \'c\' to continue...')
v.start()

# Mesh to test connection
print('Computing the mesh...')
the_shape = GroupAPI.prepare_shape_to_mesh()
the_gen = MeshGen()
the_mesh = MeshGen.create_mesh(the_gen)
the_mesh.shape_to_mesh(the_shape)

hyp2d = NetgenSimple2D(the_gen, 2.)
alg2d = NetgenAlgo2D(the_gen)
the_mesh.add_hypothesis(hyp2d, the_shape)
the_mesh.add_hypothesis(alg2d, the_shape)

the_gen.compute(the_mesh, the_shape)

v.clear()
v.display_mesh(the_mesh.object)
v.start()
