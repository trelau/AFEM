from __future__ import print_function

from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart
from asap.topology import ShapeTools

# Import model
fn = r'.\test_io\777-200LR_mod_vsp350_sref.stp'
ImportVSP.step_file(fn)

# Get components.
wing = ImportVSP.get_body('Wing')
gear = ImportVSP.get_body('Gear Pod')
fuselage = ImportVSP.get_body('Fuselage')
vtail = ImportVSP.get_body('Vtail')
htail = ImportVSP.get_body('Htail')
other_wing = ImportVSP.get_body('Wing.2')
other_htail = ImportVSP.get_body('Htail.1')

# Floor
pln = CreateGeom.plane_by_normal([0, 0, -24], [0, 0, 1])
floor = CreatePart.floor.by_sref('floor', fuselage, pln)

# Fwd bulkhead
pref = wing.eval(0.15, 0.)
pln1 = CreateGeom.plane_by_normal(pref, [1, 0, 0])
bh1 = CreatePart.floor.by_sref('fwd bulkhead', fuselage, pln1)

# Rear bulkhead
pref = wing.eval(0.65, 0.)
pln2 = CreateGeom.plane_by_normal(pref, [1, 0, 0])
bh2 = CreatePart.floor.by_sref('rear bulkhead', fuselage, pln2)

# Frames between bulkheads.
plns = CreateGeom.planes_between_planes(pln1, pln2, 24.)
frames = []
for pln in plns:
    frame = CreatePart.frame.by_sref('frame', fuselage, pln, 3.)
    frames.append(frame)

# Skin
skin = CreatePart.surface_part('fuse_skin', fuselage.shell)
skin.set_shape(fuselage.shell)

# Join
floor.join(bh1, bh2, *frames)
skin.join(floor, bh1, bh2, *frames)

fwd = ShapeTools.create_halfspace(bh1, [0, 0, 0])
aft = ShapeTools.create_halfspace(bh2, [1e6, 0, 0])
skin.discard(fwd)
skin.discard(aft)
floor.discard(fwd)
floor.discard(aft)

above_floor = ShapeTools.create_halfspace(floor, [0, 0, 1e6])

bh1.discard(above_floor)
bh2.discard(above_floor)

# View
skin.set_color(0.5, 0.5, 0.5)
skin.set_transparency(0.5)
Viewer.add_items(*[floor, skin, bh1, bh2])
Viewer.add_items(*frames)
Viewer.show()

# Mesh
for part in [floor, bh1, bh2, skin] + frames:
    part.mesh(maxh=6., quad_dominated=False)
    Viewer.add_meshes(part)
Viewer.show_mesh()
