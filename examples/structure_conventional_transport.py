from __future__ import print_function

from OCCT.BRepPrimAPI import BRepPrimAPI_MakeSphere

from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.misc.check import pairwise
from afem.oml import Body
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Set units to inch.
Settings.set_units('in')

# Import model
fn = r'../models/777-200LR.xbf'
bodies = Body.load_bodies(fn)

# Get components.
wing = bodies['Wing']
gear = bodies['Gear Pod']
fuselage = bodies['Fuselage']
vtail = bodies['Vtail']
htail = bodies['Htail']
other_wing = bodies['Wing.2']
other_htail = bodies['Htail.1']

# FUSELAGE --------------------------------------------------------------------
fuse_group = GroupAPI.create_group('fuselage group')

# Nose bulkhead
pln = PlaneByAxes((36, 0, 0), 'yz').plane
bh1 = BulkheadBySurface('bh1', pln, fuselage).bulkhead

# Flight deck bulkhead
pln = PlaneByAxes((120, 0, -24), 'yz').plane
bh2 = BulkheadBySurface('bh2', pln, fuselage).bulkhead

# Cut door in flight deck
face = FaceByPlane(pln, -18, 18, 0, 72).face
bh2.cut(face)

# Rear pressure bulkhead
p = vtail.eval(0, 0)
pln = PlaneByAxes(p, 'yz').plane
section = IntersectShapes(fuselage.shape, pln).shape
cg = LinearProps(section).cg
sphere = BRepPrimAPI_MakeSphere(cg, 120).Shell()
rev_fuselage = fuselage.shape.Reversed()
sphere = CutShapes(sphere, rev_fuselage).shape
rear_pressure_bh = Bulkhead('rear pressure bh', sphere)
hs = HalfspaceBySurface(pln, (0, 0, 0)).solid
rear_pressure_bh.discard_by_solid(hs)
rear_pressure_bh.fix()
section = IntersectShapes(rear_pressure_bh.shape, fuselage.shape).shape
cg = LinearProps(section).cg
rear_pressure_bh_pln = PlaneByAxes(cg, 'yz').plane

# Fwd bulkhead at 25% wing chord
p0 = wing.eval(0.25, 0.)
pln = PlaneByAxes(p0, 'yz').plane
fwd_bh = BulkheadBySurface('fwd bh', pln, fuselage).bulkhead

# Rear bulkhead at 65% wing chord
p0 = wing.eval(0.65, 0.)
pln = PlaneByAxes(p0, 'yz').plane
rear_bh = BulkheadBySurface('rear bh', pln, fuselage).bulkhead

# Aft bulkhead at 85% wing chord
p0 = wing.eval(0.85, 0.)
pln = PlaneByAxes(p0, 'yz').plane
aft_bh = BulkheadBySurface('aft bh', pln, fuselage).bulkhead

# Floor
pln = PlaneByAxes((0., 0., -24.), 'xy').plane
floor = FloorBySurface('floor', pln, fuselage).floor
floor.set_transparency(0.5)
floor.cut(rear_pressure_bh)
floor.discard_by_solid(hs.reversed())
floor.cut(bh1.plane)
floor.cut(bh2.plane)

# Skin
fskin = SkinByBody('fuselage skin', fuselage).skin
fskin.set_transparency(0.5)

# Cutting solids.
p0 = wing.eval(0., 0.)
fwd_pln = PlaneByAxes(p0, 'yz').plane
# fwd_cut = SolidByPlane(fwd_pln, 1.e6, 1.e6, -1.e6).solid

p0 = wing.eval(1., 0.)
aft_pln = PlaneByAxes(p0, 'yz').plane
# aft_cut = SolidByPlane(aft_pln, 1.e6, 1.e6, 1.e6).solid

shell = floor.make_shell()
above_floor = HalfspaceByShape(shell, (0., 0., 1.e6)).solid
below_floor = HalfspaceByShape(shell, (0., 0., -1.e6)).solid
above_floor2 = ShapeByDrag(shell, (0, 0, 500)).shape

left_pln = PlaneByAxes(axes='xz').plane
left_cut = SolidByPlane(left_pln, 1.e6, 1.e6, -1.e6).solid

# Frames
plns = [bh1.plane, bh2.plane, fwd_bh.sref, rear_bh.sref, aft_bh.sref,
        rear_pressure_bh_pln]
frames = []
next_index = 1
for pln1, pln2 in pairwise(plns):
    builder = FramesBetweenPlanesByDistance('frame', pln1, pln2, 24.,
                                            fuselage, 4., 24., -24.,
                                            first_index=next_index)
    frames += builder.frames
    next_index = builder.next_index

# Floor beams and posts
pln = PlaneByAxes((0, 80, 0), 'xz').plane
post_face1 = FaceByPlane(pln, -10000, 10000, -10000, 10000).face
pln = PlaneByAxes((0, -80, 0), 'xz').plane
post_face2 = FaceByPlane(pln, -10000, 10000, -10000, 10000).face
for frame in frames[3:]:
    section = IntersectShapes(floor.shape, frame.plane).shape
    shape = ShapeByDrag(section, (0, 0, -4)).shape
    shape = CutShapes(shape, rev_fuselage).shape
    beam = SurfacePart('beam', shape)
    floor.split(section, False)

    section = IntersectShapes(frame.plane, post_face1).shape
    shape = CutShapes(section, rev_fuselage).shape
    if not shape.is_null:
        post = Beam1D('post', shape)
        post.set_color(1, 0, 0)
        post.cut(floor)
        post.cut(above_floor2)

    section = IntersectShapes(frame.plane, post_face2).shape
    shape = CommonShapes(section, fuselage.shape).shape
    if not shape.is_null:
        post = Beam1D('post', shape)
        post.set_color(1, 0, 0)
        post.cut(floor)
        post.cut(above_floor2)

# Frames at bulkheads
plns = [fwd_bh.plane, rear_bh.plane, aft_bh.plane]
indx = len(frames) + 1
bh_frames = FramesByPlanes('frame', plns, fuselage, 4., indx).frames

CutParts(bh_frames, below_floor)
frames += bh_frames

# Cut and discard to get half model
# CutParts([fskin, floor], fwd_cut)
# CutParts([fskin, floor], aft_cut)
CutParts([fwd_bh, rear_bh, aft_bh], above_floor)

# Cut a piece of the wing to cut fuselage faces inside it
cut1 = SolidByPlane(fwd_bh.plane, 1.e6, 1.e6, -1.e6).solid
cut2 = SolidByPlane(aft_bh.plane, 1.e6, 1.e6, 1.e6).solid
joined_wing = FuseShapes(wing.shape, other_wing.shape).shape
bop = CutShapes()
bop.set_args([joined_wing])
bop.set_tools([cut1, cut2])
bop.build()
cutter = bop.shape
CutParts([fskin, fwd_bh, rear_bh, aft_bh, floor] + frames, cutter)

# Tail bulkhead
bbox = BBox()
bbox.add_shape(fuselage.shape)
xmax = bbox.xmax
pln = PlaneByAxes((xmax - 36, 0, 0), 'yz').plane
tail_bh = BulkheadBySurface('tail bh', pln, fuselage).bulkhead

# Aft frames
plns = [rear_pressure_bh_pln, tail_bh.plane]
frames = []
for pln1, pln2 in pairwise(plns):
    builder = FramesBetweenPlanesByDistance('frame', pln1, pln2, 24.,
                                            fuselage, 8, 36., -36.,
                                            first_index=next_index)
    frames += builder.frames
    next_index = builder.next_index

# Join
# FuseSurfaceParts([fwd_bh, rear_bh, aft_bh, floor, fskin], frames)

master = GroupAPI.get_master()
v = Viewer(1280, 1024)
v.add(master)
v.start()
