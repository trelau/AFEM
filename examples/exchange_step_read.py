"""
Test case from: https://github.com/tpaviot/pythonocc-core/issues/522
"""
from afem.exchange import StepRead
from afem.graphics import Viewer

# Read the file
reader = StepRead('../models/geometry_names.step')
shape = reader.shape

# Traverse all the faces and find the desired ones based on their name
faces = []
for f in shape.faces:
    name = reader.name_from_shape(f)
    if name == 'CYLINDER_TOP':
        f.set_color(1, 0, 0)
    elif name == 'FACE_UP':
        f.set_color(0, 0, 1)
    else:
        f.set_color(0.5, 0.5, 0.5)
    faces.append(f)

# Show the model with named faces as different colors
gui = Viewer()
gui.add(*faces)
gui.start()
