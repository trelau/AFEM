from __future__ import print_function

import time

from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.misc.check import pairwise
from afem.oml import Body
from afem.smesh import *
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Set units to inch.
Settings.set_units('in')

# Import model
fn = '../models/777-200LR.xbf'
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

# Skin
fskin = SkinByBody('fuselage skin', fuselage).skin

# Cutting solids.
p0 = wing.eval(0., 0.)
fwd_pln = PlaneByAxes(p0, 'yz').plane
fwd_cut = SolidByPlane(fwd_pln, 1.e6, 1.e6, -1.e6).solid

p0 = wing.eval(1., 0.)
aft_pln = PlaneByAxes(p0, 'yz').plane
aft_cut = SolidByPlane(aft_pln, 1.e6, 1.e6, 1.e6).solid

shell = floor.make_shell()
above_floor = HalfspaceByShape(shell, (0., 0., 1.e6)).solid
below_floor = HalfspaceByShape(shell, (0., 0., -1.e6)).solid

left_pln = PlaneByAxes(axes='xz').plane
left_cut = SolidByPlane(left_pln, 1.e6, 1.e6, -1.e6).solid

# Frames
plns = [fwd_pln, fwd_bh.sref, rear_bh.sref, aft_bh.sref, aft_pln]
frames = []
next_index = 1
for pln1, pln2 in pairwise(plns):
    builder = FramesBetweenPlanesByDistance('frame', pln1, pln2, 24.,
                                            fuselage, 4., 24., -24.,
                                            first_index=next_index)
    frames += builder.frames
    next_index = builder.next_index

# Frames at bulkheads
plns = [fwd_bh.plane, rear_bh.plane, aft_bh.plane]
indx = len(frames) + 1
bh_frames = FramesByPlanes('frame', plns, fuselage, 4., indx).frames

CutParts(bh_frames, below_floor)
frames += bh_frames

# Cut and discard to get half model
CutParts([fskin, floor], fwd_cut)
CutParts([fskin, floor], aft_cut)
CutParts([fwd_bh, rear_bh, aft_bh], above_floor)

# Cut a piece of the wing to cut fuselage faces inside it
cut1 = SolidByPlane(fwd_bh.plane, 1.e6, 1.e6, -1.e6).solid
cut2 = SolidByPlane(aft_bh.plane, 1.e6, 1.e6, 1.e6).solid
joined_wing = FuseShapes(wing.solid, other_wing.solid).shape
bop = CutShapes()
bop.set_args([joined_wing])
bop.set_tools([cut1, cut2])
bop.build()
cutter = bop.shape
CutParts([fskin, fwd_bh, rear_bh, aft_bh, floor] + frames, cutter)

# Join
FuseSurfaceParts([fwd_bh, rear_bh, aft_bh, floor, fskin], frames)

# WING ------------------------------------------------------------------------
wing_group = GroupAPI.create_group('wing group')

# Root rib

# Find the intersection between the fuselage and wing
shape = IntersectShapes(fuselage.solid, wing.solid).shape

# Discard parts of the intersection curve fwd/aft of bulkheads and keep only
# the upper skin intersection.
cut1 = SolidByPlane(fwd_bh.plane, 500., 500., -500.).solid
cut2 = SolidByPlane(rear_bh.plane, 500., 500., 500.).solid
cut3 = HalfspaceBySurface(wing.sref, (0., 0., -1000.)).solid

# Use consecutive cuts since BOPs struggle with half-spaces.
shape = CutShapes(shape, cut1).shape
shape = CutShapes(shape, cut2).shape
shape = CutShapes(shape, cut3).shape

# Make face in z-direction
rshape = ShapeByDrag(shape, (0., 0., -100.)).shape

root_rib = RibByShape('root rib', rshape, wing).rib

# Cut the root rib where the frame reference planes are. This 1) Resolves
# some issues the mesher was having with the root rib topology and 2) Allow
# for 1-D elements to be created at frame locations for stiffeners.
for frame in frames:
    root_rib.cut(frame.sref)

# Center spars between root chord and root rib
root_chord = wing.sref.v_iso(wing.sref.v1)

fc_spar = SparBetweenShapes('fc spar', root_chord, root_rib.shape, wing,
                            fwd_bh.sref).spar
rc_spar = SparBetweenShapes('rc spar', root_chord, root_rib.shape, wing,
                            rear_bh.sref).spar

# Tip rib
v = wing.sref.vknots[-2]
tip_rib = RibByParameters('tip rib', 0.15, v, 0.65, v, wing).rib

# Front and rear spar. Use intersection of root rib and center spars to
# define planes so the edges all align at the root rib.
pln = PlaneByIntersectingShapes(root_rib.shape, fc_spar.shape,
                                tip_rib.p1).plane

p1 = root_rib.p1
p2 = tip_rib.p1
fspar = SparByPoints('fspar', p1, p2, wing, pln).spar

pln = PlaneByIntersectingShapes(root_rib.shape, rc_spar.shape,
                                tip_rib.p2).plane
p1 = root_rib.p2
p2 = tip_rib.p2
rspar = SparByPoints('rspar', p1, p2, wing, pln).spar

# Aux spar
v = wing.sref.vknots[1]
p0 = wing.eval(0.5, v)
pln = PlaneByAxes(p0, 'xz').plane
aux_spar = SparBetweenShapes('aux spar', root_chord, pln, wing,
                             aft_bh.sref).spar

# Rib at end of aux spar
pln = PlaneByAxes(aux_spar.p2, 'xz').plane
face1 = FaceBySurface(fspar.sref).face
face2 = FaceBySurface(aft_bh.sref).face
kink_rib = RibBetweenShapes('kink rib', face1, face2, wing, pln).rib

# Streamwise ribs between root and kink rib.
root_pln = PlaneByAxes(root_rib.p2, 'xz').plane
plns = [root_pln, kink_rib.sref]
inbd_ribs = RibsBetweenPlanesByDistance('inbd rib', root_pln, kink_rib.plane,
                                        30., fspar.shape, rspar.shape, wing,
                                        30., -30.).ribs

# Outboard ribs
u1 = rspar.invert_cref(kink_rib.p2)
outbd_ribs = RibsAlongCurveByDistance('outbd rib', rspar.cref, 30.,
                                      fspar.shape, rspar.shape, wing, u1=u1,
                                      d1=30., d2=-30.).ribs

# Rib along kink rib to front spar.
u1 = fspar.invert_cref(kink_rib.p1)
u2 = fspar.invert_cref(outbd_ribs[0].p1)
kink_ribs = RibsAlongCurveByDistance('kink rib', fspar.cref, 30., fspar.shape,
                                     kink_rib.shape, wing, u1=u1, u2=u2,
                                     d1=30., d2=-30., nmin=1).ribs

# Center ribs.
xz_pln = PlaneByAxes(axes='xz').plane
center_ribs = RibsAlongCurveByNumber('center rib', fc_spar.cref, 3,
                                     fc_spar.shape, rc_spar.shape, wing,
                                     ref_pln=xz_pln, d1=30., d2=-30.).ribs

# Fuse wing internal structure and discard faces
internal_parts = GroupAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

# Wing skin
wskin = SkinByBody('wing skin', wing).skin

# Join the wing skin and internal structure
wskin.fuse(*internal_parts)

# Discard faces touching reference surface.
wskin.discard_by_dmin(wing.sref_shape, 1.0)

# Fix skin since it's not a single shell anymore, but a compound of two
# shells (upper and lower skin).
wskin.fix()

# JOIN ------------------------------------------------------------------------
print('Fusing assemblies...')
start = time.time()
bop = FuseGroups([wing_group, fuse_group])
print('complete', time.time() - start)

# Mesh
shape_to_mesh = GroupAPI.prepare_shape_to_mesh()
the_gen = MeshGen()
the_mesh = the_gen.create_mesh(shape_to_mesh)

# Use a single global hypothesis based on maximum length.
hyp1d = MaxLength1D(the_gen, 4.)
alg1d = Regular1D(the_gen)
the_mesh.add_hypotheses([hyp1d, alg1d], shape_to_mesh)

# Netgen unstructured quad-dominated
netgen_hyp = NetgenSimple2D(the_gen, 4.)
netgen_alg = NetgenAlgo2D(the_gen)
the_mesh.add_hypotheses([netgen_hyp, netgen_alg], shape_to_mesh)

# Apply mapped quadrangle to internal structure
mapped_hyp = QuadrangleHypo2D(the_gen)
mapped_algo = QuadrangleAlgo2D(the_gen)
for part_ in internal_parts:
    for face in part_.faces:
        if mapped_algo.is_applicable(face):
            the_mesh.add_hypotheses([mapped_algo, mapped_hyp], face)

# Compute the mesh
mesh_start = time.time()
print('Computing mesh...')
status = the_gen.compute(the_mesh, shape_to_mesh)
if not status:
    print('Failed to compute mesh')
else:
    print('Meshing complete in ', time.time() - mesh_start, ' seconds.')

# Find free edges
tool = ExploreFreeEdges(shape_to_mesh)

# View
v = Viewer()
v.display_mesh(the_mesh.object, 2)
v.add(*tool.free_edges)
v.start()
