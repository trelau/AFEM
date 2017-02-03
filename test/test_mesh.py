from __future__ import print_function

import time

from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart, PartTools

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

start = time.time()
# Spars
fspar = CreatePart.spar.by_parameters('fspar', wing, 0.15, 0.05, 0.15, 0.925)
rspar = CreatePart.spar.by_parameters('rspar', wing, 0.65, 0.05, 0.65, 0.925)

# Root and tip rib.
root = CreatePart.rib.by_points('root rib', wing, fspar.p1, rspar.p1)
tip = CreatePart.rib.by_points('tip rib', wing, fspar.p2, rspar.p2)
ribs = [root]

# Ribs
p2 = rspar.distribute_points(30., 30., -30.)
p1 = [p.copy() for p in p2]
fspar.project_points(p1)
for pf, pr in zip(p1, p2):
    rib = CreatePart.rib.by_points('rib', wing, pf, pr)
    if not rib:
        continue
    ribs.append(rib)
ribs.append(tip)

# Fuse internal structure
status = PartTools.join_wing_parts([fspar, rspar] + ribs)

# Discard faces on internal structure
for part in [fspar, rspar] + ribs:
    part.discard()

# Skin
skin = CreatePart.surface_part('wing skin', wing.shell)
skin.set_shape(wing.shell)

# Join the fuse_skin and internal structure
skin.join(fspar, rspar, *ribs)
skin.fix()

# Discard faces touching reference surface.
skin.discard(wing.sref, 0.005)

# Test a mesh.
for part in [fspar, rspar, skin] + ribs:
    part.mesh(6., quad_dominated=False)
    Viewer.add_meshes(part)
Viewer.show_mesh()
