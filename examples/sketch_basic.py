from afem.geometry import *
from afem.graphics import Viewer
from afem.sketch import *
from afem.topology import *

# Create a new cross section
cs = Airfoil()

# Generate a 2-D profile by reading and approximating an airfoil file from the
# UIUC database. Close the trailing edge if necessary.
cs.read_uiuc('../models/clarky.dat', close=True)

# Define a plane at the root and scale
pln1 = PlaneByAxes(axes='xz').plane
cs.build(pln1, scale=5)
wire1 = cs.wires[0]

# Define plane at the tip and rotate
pln2 = PlaneByAxes((3, 15, 0), axes='xz').plane
cs.build(pln2, scale=1.5, rotate=3)
wire2 = cs.wires[0]

# Use the wires to loft a solid
shape = LoftShape([wire1, wire2], True).shape

gui = Viewer()
gui.add(wire1, wire2, shape)
gui.start()
