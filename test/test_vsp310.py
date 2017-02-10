from __future__ import print_function

from OCC.ShapeAnalysis import ShapeAnalysis_ShapeTolerance

from asap.io import ImportVSP
from asap.graphics import Viewer

fn = r'.\test_io\777-200LR_mod_vsp310_split.stp'

ImportVSP.step_file(fn)

bodies = ImportVSP.get_bodies()
for name in bodies:
    b = bodies[name]
    print('Body name/type ID/tolerance:', name, b.ShapeType(),
          ShapeAnalysis_ShapeTolerance().Tolerance(b, 0))
    Viewer.add_items(b)
Viewer.show()

body1 = ImportVSP.get_body('Body.1')
body2 = ImportVSP.get_body('Body.2')

body = body1.fuse(body2)
print(body)
Viewer.add_items(body)
Viewer.show()
