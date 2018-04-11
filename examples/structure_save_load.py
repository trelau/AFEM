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
wing = vsp_import.get_body('WingGeom')

# Build structure
wingbox = GroupAPI.create_group('wing box')
fspar = SparByParameters('front spar', 0.15, 0., 0.15, 1., wing).part
rspar = SparByParameters('rear spar', 0.70, 0., 0.70, 1., wing).part
RibByPoints('root rib', fspar.p1, rspar.p1, wing)
RibByPoints('tip rib', fspar.p2, rspar.p2, wing)
RibsAlongCurveByDistance('rib', rspar.cref, 30, fspar.shape, rspar.shape,
                         wing, d1=30, d2=-30)
internal_parts = wingbox.get_parts()
skin = SkinByBody('skin', wing).part
cref = wing.sref.u_iso(0.5)
skin.discard_by_dmin(cref, 1.0)

FuseSurfaceParts([skin], internal_parts)

group = GroupAPI.get_master()

v = Viewer()
v.add(group)
v.start()

# Save
print('\nSaving group...')
GroupAPI.save_model('structure.xbf')
print('done.\n')

# Load into new group
print('Loading group...')
new_group = GroupAPI.create_group('new model')
GroupAPI.load_model('structure.xbf', new_group)
print('done.')
v.clear()
v.add(new_group)
v.start()
