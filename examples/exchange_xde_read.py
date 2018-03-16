from afem.exchange.xde import XdeRead
from afem.graphics import Viewer

fn = '../models/777-200LR_Onshape.step'

reader = XdeRead(fn)

v = Viewer()
for label in reader.labels:
    shape = reader.get_shape(label)
    print('{}: {}'.format(label, shape))
    v.add(shape)
v.start()
