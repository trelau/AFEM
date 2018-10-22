from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer
from afem.structure import *

Settings.log_to_console()

# Set units to inch.
Settings.set_units('in')

# Import model
fn = r'../models/simple_wing.stp'
vsp_import = ImportVSP(fn)
wing = vsp_import['WingGeom']

# Build structure
wingbox = GroupAPI.create_group('wing box')
fspar = SparByParameters('front spar', 0.15, 0., 0.15, 1., wing).part
rspar = SparByParameters('rear spar', 0.70, 0., 0.70, 1., wing).part
RibByPoints('root rib', fspar.cref.p1, rspar.cref.p1, wing)
RibByPoints('tip rib', fspar.cref.p2, rspar.cref.p2, wing)
RibsAlongCurveByDistance('rib', rspar.cref, 30, fspar.shape, rspar.shape,
                         wing, d1=30, d2=-30)
internal_parts = wingbox.get_parts()
skin = SkinByBody('skin', wing).part
cref = wing.sref.u_iso(0.5)
skin.discard_by_dmin(cref, 1.0)

FuseSurfaceParts([skin], internal_parts)

# Mesh
print('Creating mesh')
mesh = MeshVehicle(4.)
mesh.compute()

# View
gui = Viewer()
gui.add(wingbox)
gui.start()
gui.clear()
gui.add(mesh.mesh)
gui.start()
