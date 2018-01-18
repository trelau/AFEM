"""
References:

    * https://github.com/ramcdona/Code-Eli/commit/4fb2abc9a67d5b8e335da16e57b09b3db19df8d5

    * https://github.com/OpenVSP/OpenVSP/commit/a48c1daf0616cecbd0efef456f5601c075f70adc

"""

from afem.config import Settings
from afem.exchange import StepRead
from afem.geometry import *
from afem.graphics import Viewer
from afem.topology import *

# Set units to inch
Settings.set_units('ft')

# OpenVSP 3.5.0 ----------------------------------------------------------------

fn = r'../models/DefaultWing-OpenVSP-3-5-0.stp'

step = StepRead(fn)
shape1 = step.shape

face1 = ExploreShape.get_faces(shape1)[0]
surf1 = ExploreShape.surface_of_face(face1)
surf1 = NurbsSurface.downcast(surf1)
print('---OpenVSP 3.5.0 Default Wing Surface---')
print('\tp={}, q={}, nuknots={}, nvknots={}, area={} sq. ft'.format(surf1.p, surf1.q, len(surf1.uknots), len(surf1.vknots), surf1.area))

v = Viewer()
v.add(shape1)
v.start()

# OpenVSP 3.15.0 ---------------------------------------------------------------

fn = r'../models/DefaultWing-OpenVSP-3-15-0.stp'

step = StepRead(fn)
shape2 = step.shape

face2 = ExploreShape.get_faces(shape2)[0]
surf2 = ExploreShape.surface_of_face(face2)
surf2 = NurbsSurface.downcast(surf2)
print('\n---OpenVSP 3.15.0 Default Wing Surface---')
print('\tp={}, q={}, nuknots={}, nvknots={}, area={} sq. ft'.format(surf2.p, surf2.q, len(surf2.uknots), len(surf2.vknots), surf2.area))

v.add(shape2)
v.start()

# Split the shapes to make C1 continuous ---------------------------------------
print('\n---Splitting OpenVSP 3.5.0 shape to make C1 continuous---')
v.clear()
c0_shape1 = DivideC0Shape(shape1).shape
c0_shape2 = DivideC0Shape(shape2).shape
v.add(c0_shape1)
v.start()
v.clear()
print('\n---Splitting OpenVSP 3.15.0 shape to make C1 continuous---')
v.add(c0_shape2)
v.start()
