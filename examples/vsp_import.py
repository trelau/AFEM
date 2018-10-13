from afem.config import Settings
from afem.exchange import ImportVSP
from afem.graphics import Viewer

Settings.log_to_console()

fn = '../models/777-200LR.stp'

vsp_import = ImportVSP(fn)

gui = Viewer()
gui.add(*vsp_import.all_bodies)
gui.start()
