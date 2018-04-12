from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer

Settings.log_to_console()

fn = '../models/777-200LR.stp'

vsp_import = ImportVSP(fn)

v = Viewer()
v.add(*vsp_import.all_bodies)
v.start()
