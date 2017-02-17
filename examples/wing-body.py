from __future__ import print_function

import time

from asap.geometry import CreateGeom
from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import AssemblyMgr, CreatePart, PartTools
from asap.topology import ShapeTools

# Import model
fn = './models/777-200LR_mod_vsp350_sref.stp'
ImportVSP.step_file(fn)

# Get components.
wing = ImportVSP.get_body('Wing')
gear = ImportVSP.get_body('Gear Pod')
fuselage = ImportVSP.get_body('Fuselage')
vtail = ImportVSP.get_body('Vtail')
htail = ImportVSP.get_body('Htail')
other_wing = ImportVSP.get_body('Wing.2')
other_htail = ImportVSP.get_body('Htail.1')

# FUSELAGE --------------------------------------------------------------------
AssemblyMgr.create_assy('fuselage assy')

# Fwd bulkhead at 25% wing chord
p0 = wing.eval(0.25, 0.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
fwd_bh = CreatePart.bulkhead.by_sref('fwd bh', fuselage, pln)

# Rear bulkhead at 65% wing chord
p0 = wing.eval(0.65, 0.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
rear_bh = CreatePart.bulkhead.by_sref('rear bh', fuselage, pln)

# Aft bulkhead at 85% wing chord
p0 = wing.eval(0.85, 0.)
pln = CreateGeom.plane_by_axes(p0, 'yz')
aft_bh = CreatePart.bulkhead.by_sref('aft bh', fuselage, pln)

# Floor
pln = CreateGeom.plane_by_axes([0, 0, -24], 'xy')
floor = CreatePart.floor.by_sref('floor', fuselage, pln)

# Skin
fskin = CreatePart.skin.from_body('fuselage skin', fuselage)

# Cutting solids.
p0 = wing.eval(0., 0.)
fwd_pln = CreateGeom.plane_by_axes(p0, 'yz')
fwd_cut = ShapeTools.box_from_plane(fwd_pln, 1e6, 1e6, -1e6)

p0 = wing.eval(1., 0.)
aft_pln = CreateGeom.plane_by_axes(p0, 'yz')
aft_cut = ShapeTools.box_from_plane(aft_pln, 1e6, 1e6, 1e6)

above_floor = ShapeTools.make_halfspace(floor, [0, 0, 1e6])
below_floor = ShapeTools.make_halfspace(floor, [0, 0, -1e6])

left_pln = CreateGeom.plane_by_axes([0., 0., 0.], 'xz')
left_cut = ShapeTools.box_from_plane(left_pln, 1e6, 1e6, -1e6)

# Frames
plns = [fwd_pln, fwd_bh.sref, rear_bh.sref, aft_bh.sref, aft_pln]
frames = CreatePart.frame.between_planes('frame', fuselage, plns, 4., 24.)

# Frames at bulkheads
plns = [fwd_bh.sref, rear_bh.sref, aft_bh.sref]
indx = len(frames) + 1
bh_frames = CreatePart.frame.at_shapes('frame', fuselage, plns, 4., indx)
PartTools.cut_parts(bh_frames, below_floor)
frames += bh_frames

# Cut and discard to get half model
PartTools.cut_parts([fskin, floor], fwd_cut)
PartTools.cut_parts([fskin, floor], aft_cut)
PartTools.cut_parts([fwd_bh, rear_bh, aft_bh], above_floor)

# Cut a piece of the wing to cut fuselage faces inside it.
cut1 = ShapeTools.box_from_plane(fwd_bh.sref, 1e6, 1e6, -1e6)
cut2 = ShapeTools.box_from_plane(aft_bh.sref, 1e6, 1e6, 1e6)
joined_wing = wing.fuse(other_wing)
cutter = joined_wing.bop_algo([cut1, cut2], 'cut')
PartTools.cut_parts([fskin, fwd_bh, rear_bh, aft_bh, floor] + frames, cutter)

# Join
PartTools.fuse_parts([fwd_bh, rear_bh, aft_bh, floor, fskin] + frames)

# WING ------------------------------------------------------------------------
AssemblyMgr.create_assy('wing assy')

# Root rib

# Find the intersection between the fuselage and wing
shape = fuselage.section(wing)

# Discard parts of the intersection curve fwd/aft of bulkheads and keep only
# the upper skin intersection.
cut1 = ShapeTools.box_from_plane(fwd_bh.sref, 500, 500, -500)
cut2 = ShapeTools.box_from_plane(rear_bh.sref, 500, 500, 500)
cut3 = ShapeTools.make_halfspace(wing.sref, [0, 0, -1000])

# Use consecutive cuts since BOPs struggle with half-spaces.
shape = ShapeTools.bcut(shape, cut1)
shape = ShapeTools.bcut(shape, cut2)
shape = ShapeTools.bcut(shape, cut3)

# Make face in z-direction
vz = CreateGeom.vector([0, 0, -100])
rshape = ShapeTools.make_prism(shape, vz)

root_rib = CreatePart.rib.by_sref('root rib', wing, rshape)

# Center spars between root chord and root rib.
root_chord = wing.isocurve(v=wing.v1)

fc_spar = CreatePart.spar.between_geom('fc spar', wing, root_chord,
                                       root_rib.cref, fwd_bh.sref)

rc_spar = CreatePart.spar.between_geom('rc spar', wing, root_chord,
                                       root_rib.cref, rear_bh.sref)

# Tip rib
v = wing.vknots[-2]
tip_rib = CreatePart.rib.by_parameters('tip rib', wing, 0.15, v, 0.65, v)

# Front and rear spar. Use intersection of root rib and center spars to
# define planes so the edges all align at the root rib.
pln = ShapeTools.plane_from_section(root_rib.rshape, fc_spar.rshape,
                                    tip_rib.p1)
p1 = root_rib.p1
p2 = tip_rib.p1
fspar = CreatePart.spar.by_points('fspar', wing, p1, p2, pln)

pln = ShapeTools.plane_from_section(root_rib.rshape, rc_spar.rshape,
                                    tip_rib.p2)
p1 = root_rib.p2
p2 = tip_rib.p2
rspar = CreatePart.spar.by_points('rspar', wing, p1, p2, pln)

# Aux spar
v = wing.vknots[1]
p0 = wing.eval(0.5, v)
pln = CreateGeom.plane_by_axes(p0, 'xz')
aux_spar = CreatePart.spar.between_geom('aux spar', wing, root_chord,
                                        pln, aft_bh.sref)

# Rib at end of aux spar
pln = CreateGeom.plane_by_axes(aux_spar.p2, 'xz')
kink_rib = CreatePart.rib.between_geom('kink rib', wing, fspar.sref,
                                       aft_bh.sref, pln)

# Streamwise ribs between root and kink rib.
root_pln = CreateGeom.plane_by_axes(root_rib.p2, 'xz')
plns = [root_pln, kink_rib.sref]
inbd_ribs = CreatePart.rib.between_planes('rib', wing, plns, fspar.sref,
                                          rspar.sref, 30.)

# Outboard ribs
u1 = rspar.invert(kink_rib.p2)
outbd_ribs = CreatePart.rib.along_curve('outbd rib', wing, rspar.cref,
                                        fspar.sref, rspar.sref, 30., u1=u1,
                                        s1=30., s2=-30.)

# Rib along kink rib to front spar.
u1 = fspar.invert(kink_rib.p1)
u2 = fspar.invert(outbd_ribs[0].p1)
kink_ribs = CreatePart.rib.along_curve('kink rib', wing, fspar.cref,
                                       fspar.sref, kink_rib.sref, maxd=30.,
                                       npts=1, u1=u1, u2=u2, s1=30., s2=-30.)

# Center ribs.
xz_pln = CreateGeom.plane_by_axes(axes='xz')
center_ribs = CreatePart.rib.along_curve('center rib', wing, fc_spar.cref,
                                         fc_spar.sref, rc_spar.sref, npts=3,
                                         ref_pln=xz_pln, s1=30., s2=-30.)

# Fuse wing internal structure and discard faces
internal_parts = AssemblyMgr.get_parts()
PartTools.fuse_wing_parts(internal_parts)
PartTools.discard_faces(internal_parts)

# Wing skin
wskin = CreatePart.skin.from_body('wing skin', wing)

# Join the wing skin and internal structure
wskin.fuse(*internal_parts)

# Discard faces touching reference surface.
wskin.discard(wing.sref)

# Fix skin since it's not a single shell anymore, but a compound of two
# shells (upper and lower skin).
wskin.fix()

print('start cutting')
frame_cut = ShapeTools.make_compound(frames)
start = time.time()
wskin.cut(frame_cut)
PartTools.cut_parts([wskin, fc_spar, aux_spar], fuselage.shell)
print('done cutting', time.time() - start)

wing_parts = AssemblyMgr.get_parts('wing assy')
fuselage_parts = AssemblyMgr.get_parts('fuselage assy')
print('starting assy join')
start = time.time()
PartTools.sew_parts(wing_parts + fuselage_parts)
print('complete', time.time() - start)

# Viewing
wskin.set_transparency(0.5)
wskin.set_color(0.5, 0.5, 0.5)

for part in wing_parts + fuselage_parts:
    part.mesh(maxh=4., quad_dominated=False)
    Viewer.add_meshes(part)
Viewer.show_mesh()

for oml in [wing, gear, fuselage]:
    oml.set_color(0.5, 0.5, 0.5)
    oml.set_transparency(0.5)

# Viewer.add_items(wing)
Viewer.add_items(fwd_bh, rear_bh, aft_bh, fskin, floor, *frames)
Viewer.add_items(fc_spar, rc_spar, root_rib, tip_rib, fspar, rspar)
Viewer.add_items(wskin)
Viewer.show()

# from asap.io import StepExport
#
# assy = ShapeTools.make_compound([wing_compound, fuselage_compound])
# step_out = StepExport()
# step_out.transfer(assy)
# step_out.write('wing-body.stp')
