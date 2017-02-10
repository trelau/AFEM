from __future__ import print_function

import time

from asap.geometry import CreateGeom, ProjectGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart, PartTools
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

# wingbox ---------------------------------------------------------------------
# Spars
fspar = CreatePart.spar.by_parameters('fspar', wing, 0.15, 0.05, 0.15, 0.925)
rspar = CreatePart.spar.by_parameters('rspar', wing, 0.65, 0.05, 0.65, 0.925)

# Root and tip rib.
root = CreatePart.rib.by_points('root rib', wing, fspar.p1, rspar.p1)
tip = CreatePart.rib.by_points('tip rib', wing, fspar.p2, rspar.p2)
ribs = [root]

# Ribs
p2 = rspar.distribute_points(30.)
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
wing_skin = CreatePart.surface_part('wing wing_skin', wing.shell)
wing_skin.set_shape(wing.shell)

# Fuse internal structure
wing_parts = [fspar, rspar, fc_spar, rc_spar, mspar] + ribs
status = PartTools.fuse_wing_parts(wing_parts)

# Discard faces on internal structure
for part in [fspar, rspar, fc_spar, rc_spar, mspar] + ribs:
    part.discard()

# Join the fuse_skin and internal structure
wing_skin.fuse(fspar, rspar, fc_spar, rc_spar, mspar, *ribs)
wing_skin.fix()

# Discard faces touching reference surface.
wing_skin.discard(wing.sref, 0.005)

wing_parts = ShapeTools.make_compound([wing_skin, fspar, rspar, fc_spar,
                                       rc_spar] + ribs)

# fuselage --------------------------------------------------------------------
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
    frame = CreatePart.frame.by_sref('frame', fuselage, pln, 12)
    frames.append(frame)

# Skin
fuse_skin = CreatePart.surface_part('fuse wing_skin', fuselage.shell)
fuse_skin.set_shape(fuselage.shell)

# Join
floor.fuse(bh1, bh2, *frames)
fuse_skin.fuse(floor, bh1, bh2, *frames)

fwd = ShapeTools.make_halfspace(bh1, [0, 0, 0])
aft = ShapeTools.make_halfspace(bh2, [1e6, 0, 0])
fuse_skin.discard(fwd)
fuse_skin.discard(aft)
floor.discard(fwd)
floor.discard(aft)

above_floor = ShapeTools.make_halfspace(floor, [0, 0, 1e6])

bh1.discard(above_floor)
bh2.discard(above_floor)

# print('cutting fuse parts with wing...')
# start = time.time()
# for part in [fuse_skin, bh1, bh2] + frames:
#     part.cut(wing)
# print('Complete in ', time.time() - start, ' seconds.')

fuse_parts = ShapeTools.make_compound([fuse_skin, floor, bh1, bh2] + frames)

# display.DisplayShape(wing_parts)
# display.DisplayShape(fuse_parts, update=True)
# start_display()

# Join ------------------------------------------------------------------------
print('starting wing-body fuse...')
start = time.time()
shape = ShapeTools.bfuse(wing_parts, fuse_parts)
print('Complete in ', time.time() - start, ' seconds.')

faces = ShapeTools.get_faces(shape)
for f in faces:
    Viewer.display(f, 'random', 0.5)
Viewer.show()
