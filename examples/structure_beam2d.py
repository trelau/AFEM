from afem.config import Settings
from afem.exchange import ImportVSP
from afem.geometry import *
from afem.graphics import Viewer
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Import model
fname = '../models/uniform_wing.stp'
vsp = ImportVSP(fname)
wing = vsp['Wing']

# Reference curves using wing reference surface generated from modified
# OpenVSP
cref1 = wing.sref.u_iso(0.25)
cref2 = wing.sref.u_iso(0.65)

# Plane at rear spar to cut ribs with. The will enable the 1-D beam to be
# properly fused.
aft_pln = wing.extract_plane(0.65, 0., 0.65, 1.)

# Ribs
xz_pln = PlaneByAxes(axes='xz').plane
ribs = RibsAlongCurveByDistance('Rib', cref1, 16., cref1, cref2, wing, xz_pln,
                                d1=6).parts

# Cut the ribs so there is an edge to join the 1-D beam with
CutParts(ribs, aft_pln)

# Circular cross section for front spar
circle1 = CircleByNormal(cref1.eval(0.), cref1.deriv(0., 1), 2.75).circle
beam1 = Beam2DBySweep('beam 1', cref1, circle1).part

# Aft spar is a 1-D beam
beam2 = Beam1DByCurve('beam 2', cref2).part

# Create a solid to cut the ribs with
face = FaceByPlanarWire(circle1).face
solid = SweepShape(cref1, face).shape
CutParts(ribs, solid)

# Wing skin
skin = SkinByBody('skin', wing).part
skin.set_color(0.5, 0.5, 0.5)
skin.set_transparency(0.5)

# Discard end caps of skin
skin.discard_by_dmin(cref1, 0.5)

# Join all the parts
parts = GroupAPI.get_parts()
SplitParts(parts)

gui = Viewer()
gui.add(*parts)
gui.start()

# Mesh
print('Computing the mesh...')
mesh = MeshVehicle(1.)
mesh.compute()
gui.clear()
gui.add(mesh)
gui.start()
