from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer

Settings.log_to_console()

# This model has bad shapes. Look at the leading edge and you'll see the
# self-intersection.
fn = '../models/vsp_bad_geom.stp'

vsp_import = ImportVSP(fn)
htail = vsp_import.get_body('Htail')

v = Viewer()
v.display(htail.solid, (0.5, 0.5, 0.5))
for shape in vsp_import.invalid_shapes:
    v.display(shape, (1., 0., 0.))
v.start()
