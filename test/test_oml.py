from __future__ import print_function

from OCC.ShapeAnalysis import ShapeAnalysis_ShapeTolerance

from asap.graphics import Viewer
from asap.io import ImportVSP

# fn = r'.\test_io\777-200LR_mod_vsp350_split_sref.stp'
# fn = r'.\test_io\777-200LR_mod_vsp350_sref.stp'
# fn = r'.\test_io\777-200LR_mod_vsp310_split.stp'
# fn = r'.\test_io\TBW_SUGAR_mod_vsp350_sref.stp'
# fn = r'.\test_io\F-16.stp'
# fn = r'.\test_io\HWB.stp'
# fn = r'.\test_io\N2A.stp'
fn = r'.\test_io\manta2.stp'

ImportVSP.step_file(fn)

bodies = ImportVSP.get_bodies()
for name in bodies:
    b = bodies[name]
    print('Body name/type ID/tolerance:', name, b.ShapeType(),
          ShapeAnalysis_ShapeTolerance().Tolerance(b, 0))
    Viewer.add_items(b)
Viewer.show()
