from afem.config import Settings
from afem.exchange import StepImport
from afem.graphics import Viewer
from afem.topology import *

Settings.log_to_console()

# Set units to inch
Settings.set_units('in')

# OpenVSP 3.5.0 ----------------------------------------------------------------

fn = r'../models/777-200LR_nosplit.stp'
step = StepImport()
step.read(fn)

v = Viewer()
v.set_white_background()
v.display_shape(step.shape, (0.5, 0.5, 0.5))
v.start()
v.remove_all()

c0_shape = DivideC0Shape(step.shape).shape

for f in ExploreShape.get_faces(c0_shape):
    v.display_shape(f)
v.start()
v.remove_all()

# OpenVSP 3.15.0 ---------------------------------------------------------------

fn = r'../models/777-200LR_nosplit_vsp315.stp'
step = StepImport()
step.read(fn)

v.display_shape(step.shape, (0.5, 0.5, 0.5))
v.start()
v.remove_all()

c0_shape = DivideC0Shape(step.shape).shape

for f in ExploreShape.get_faces(c0_shape):
    v.display_shape(f)
v.start()
