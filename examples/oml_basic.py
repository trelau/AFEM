from math import tan, radians

from afem.geometry import *
from afem.graphics import Viewer
from afem.oml import *
from afem.sketch import *
from afem.topology import *

# Parameters
semispan = 107.  # wing semi-span
sweep = 34.  # Leading edge sweep
uk = 0.28  # Percent semi-span to locate section cross section
c1 = 51.5  # Root chord
c2 = 31.  # Chord of second section
c3 = 7.  # Tip chord
t3 = 3.  # Tip washout in degrees

# Define leading edge points of cross sections
p1x, p1y = 0., 0.
p2x = semispan * uk * tan(radians(sweep))
p2y = semispan * uk
p3x = semispan * tan(radians(sweep))
p3y = semispan

# Create a cross section using an UIUC airfoil file
cs = Airfoil()
cs.read_uiuc('../models/clarky.dat')

# Define cross section planes
pln1 = PlaneByAxes((p1x, p1y, 0), axes='xz').plane
pln2 = PlaneByAxes((p2x, p2y, 0), axes='xz').plane
pln3 = PlaneByAxes((p3x, p3y, 0), axes='xz').plane

# Build cross sections
cs.build(pln1, scale=c1)
wire1 = cs.wires[0]
cs.build(pln2, scale=c2)
wire2 = cs.wires[0]
cs.build(pln3, scale=c3, rotate=t3)
wire3 = cs.wires[0]

# Loft a solid
shape = LoftShape([wire1, wire2, wire3], True, make_ruled=True).shape

# Make a body
wing = Body(shape, 'Wing')
wing.set_transparency(0.5)
wing.set_color(1, 0, 0)

# Build a reference surface using chord lines of the cross sections. Make sure
# to use the same scaling and rotation parameters. Set the parametric domains
# to be between 0 and 1 for convenience.
chord1 = cs.build_chord(pln1, scale=c1)
chord2 = cs.build_chord(pln2, scale=c2)
chord3 = cs.build_chord(pln3, scale=c3, rotate=t3)
sref = NurbsSurfaceByInterp([chord1, chord2, chord3], 1).surface
sref.set_udomain(0., 1.)
sref.set_vdomain(0., 1.)

# Set the wing reference surface
wing.set_sref(sref)

# Show wing and reference surface
gui = Viewer()
gui.add(wing, wing.sref)
gui.start()

# Show the underlying shape of the reference surface
gui.clear()
gui.add(wing.sref_shape)
gui.start()

# Evaluate point
p = wing.sref.eval(0.5, 0.5)

# Extract a plane
pln = wing.extract_plane(0.5, 0., 0.5, 0.5)
face = FaceByPlane(pln, -10, 10, -10, 10).face

# Extract a trimmed curve
crv = wing.extract_curve(0.15, 0.05, 0.15, 0.95)

gui.clear()
gui.add(wing.sref, p, face, crv)
gui.start()
