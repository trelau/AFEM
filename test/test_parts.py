from __future__ import print_function

from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart

fn = r'.\test_io\777-200LR_mod_vsp350_sref.stp'

ImportVSP.step_file(fn)

wing = ImportVSP.get_body('Wing')
gear = ImportVSP.get_body('Gear Pod')
fuselage = ImportVSP.get_body('Fuselage')
vtail = ImportVSP.get_body('Vtail')
htail = ImportVSP.get_body('Htail')
other_wing = ImportVSP.get_body('Wing.2')
other_htail = ImportVSP.get_body('Htail.1')

# Frame
pln = CreateGeom.plane_by_normal([1200, 0, 0], [1, 0, 0])
frame = CreatePart.frame.by_sref('frame', fuselage, pln, 2.5)

# Bulkhead
pln = CreateGeom.plane_by_normal([1000, 0, 0], [1, 0, 0])
bulkhead = CreatePart.bulkhead.by_sref('bulkhead', fuselage, pln, build=False)
bulkhead.form(wing)
bulkhead.form(other_wing)
bulkhead.build(True)

# Floor
pln = CreateGeom.plane_by_normal([1000, 0, -12], [0, 0, 1])
floor = CreatePart.floor.by_sref('floor', fuselage, pln)

# Join parts.
floor.join(bulkhead, frame)

Viewer.add_items(frame, bulkhead, floor)
Viewer.show()
