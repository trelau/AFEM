from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer

Settings.log_to_console()

# v = Viewer()
fn = r'..\models\777-200LR.stp'
# fn = r'..\models\TBW_SUGAR.stp'
# fn = r'..\models\F-16.stp'
# fn = r'..\models\HWB.stp'
# fn = r'..\models\N2A_nosplit.stp'
# fn = r'..\models\manta2.stp'
# fn = r'..\models\777-200LR_vsp315.stp'
# fn = r'..\models\777-200LR_nosplit_vsp315.stp'

ImportVSP.step_file(fn)

bodies = ImportVSP.get_bodies()
shapes = []
for name in bodies:
    b = bodies[name]
    shapes.append(b)

v = Viewer()
v.add(*shapes)
v.start()
