from afem.exchange import *
from afem.graphics import Viewer

fn = '../models/777-200LR_Onshape.step'

doc = XdeDocument()
shapes_label = doc.read_step(fn)

gui = Viewer()
for child_label in shapes_label.children_iter:
    gui.add(child_label.shape)
    print('Name: {}'.format(child_label.name))
gui.start()
