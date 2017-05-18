from asap.graphics import Viewer
from asap.io import ImportVSP, StepExport
from asap.topology import ShapeTools

fn = r'..\models\777-200LR.stp'
# fn = r'..\models\TBW_SUGAR.stp'
# fn = r'..\models\F-16.stp'
# fn = r'..\models\HWB.stp'
# fn = r'..\models\N2A_nosplit.stp'
# fn = r'..\models\manta2.stp'

ImportVSP.step_file(fn)

bodies = ImportVSP.get_bodies()
for name in bodies:
    b = bodies[name]
    print('\nBody name:', name)
    print('    Average tolerance:', ShapeTools.get_tolerance(b))
    print('    Maximum tolerance:', ShapeTools.get_tolerance(b, 1))
    Viewer.add_items(b)
print()
Viewer.show()

# Export translated VSP model to a STEP file.
compound = ShapeTools.make_compound(bodies.values())
step = StepExport(units='inch')
step.transfer(compound)
step.write('vsp2asap.stp')
