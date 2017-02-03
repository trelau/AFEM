from __future__ import print_function

from OCC.ShapeAnalysis import ShapeAnalysis_ShapeTolerance

from asap.io import ImportVSP
from asap.graphics import Viewer

# fn = r'.\test_io\777-200LR_mod_vsp350_split_sref.stp'
fn = r'.\test_io\777-200LR_mod_vsp350_sref.stp'
# fn = r'.\test_io\TBW_SUGAR_mod_vsp350_split.stp'
# fn = r'.\test_io\TBW_SUGAR_mod_vsp350_sref.stp'

ImportVSP.step_file(fn)

bodies = ImportVSP.get_bodies()
for name in bodies:
    b = bodies[name]
    print('Body name/type ID/tolerance:', name, b.ShapeType(),
          ShapeAnalysis_ShapeTolerance().Tolerance(b, 0))
    Viewer.add_items(b)
Viewer.show()

wing = ImportVSP.get_body('Wing')
gear = ImportVSP.get_body('Gear Pod')
fuselage = ImportVSP.get_body('Fuselage')
vtail = ImportVSP.get_body('Vtail')
htail = ImportVSP.get_body('Htail')
other_wing = ImportVSP.get_body('Wing.2')
other_htail = ImportVSP.get_body('Htail.1')
# other_gear = ImportVSP.get_body('Gear Pod.1')
# jury = ImportVSP.get_body('Jury')
# strut = ImportVSP.get_body('Strut')

body = fuselage.fuse(gear)
body = body.fuse(wing)
body = body.fuse(other_wing)
body = body.fuse(vtail)
body = body.fuse(htail)
body = body.fuse(other_htail)
body.set_transparency(0.5)
Viewer.add_items(body)
Viewer.show()
