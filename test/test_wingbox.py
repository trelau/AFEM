from __future__ import print_function

import time

from asap.geometry import CreateGeom, ProjectGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart, PartTools, AssemblyMgr
from asap.topology import ShapeTools

# Import model
fn = r'.\test_io\777-200LR_mod_vsp350_split_sref.stp'
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

umid = root.local_to_global(0.5)
p1 = root.eval(umid)
umid = ribs[10].local_to_global(0.5)
p2 = ribs[10].eval(umid)
mspar = CreatePart.spar.by_points('mid spar', wing, p1, p2)

# Construction geom for center structure.
root_chord = wing.extract_curve((0, 0), (1, 0))

# Front center spar.
p2 = fspar.p1
p1 = p2.copy()
ProjectGeom.point_to_geom(p1, root_chord, True)
fc_spar = CreatePart.spar.by_points('fc spar', wing, p1, p2)

# Rear center spar.
p2 = rspar.p1
p1 = p2.copy()
ProjectGeom.point_to_geom(p1, root_chord, True)
rc_spar = CreatePart.spar.by_points('rc spar', wing, p1, p2)

# Center rib
p1 = fc_spar.p1
pref = CreateGeom.plane_by_normal(p1, [0, 1, 0])
p2 = rc_spar.p1
crib = CreatePart.rib.by_points('center rib', wing, p1, p2, pref)

# Skin
skin = CreatePart.surface_part('wing skin', wing.shell)
skin.set_shape(wing.shell)

# Fuse internal structure
wing_parts = [fspar, rspar, fc_spar, rc_spar, mspar] + ribs
status = PartTools.fuse_wing_parts(wing_parts)

# Discard faces on internal structure
for part in [fspar, rspar, fc_spar, rc_spar, mspar] + ribs:
    part.discard()

# Join the fuse_skin and internal structure
skin.fuse(fspar, rspar, fc_spar, rc_spar, mspar, *ribs)
skin.fix()

print('wing_skin tol:', skin.tol)

# Discard faces touching reference surface.
skin.discard(wing.sref, 0.005)

print(time.time() - start)

# View
parts = AssemblyMgr.get_parts()
skin.set_color(0.5, 0.5, 0.5)
skin.set_transparency(0.5)
Viewer.add_items(*parts)
Viewer.show()

compound = ShapeTools.make_compound([fspar, rspar, fc_spar, rc_spar,
                                     skin, mspar] + ribs)

faces = ShapeTools.get_faces(compound)
for f in faces:
    Viewer.add_entity(f, 'random')
Viewer.show()

for rib in ribs:
    edges = rib.shared_edges(skin)
    for e in edges:
        Viewer.add_entity(e)
Viewer.show()
