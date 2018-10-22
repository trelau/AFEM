from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.oml import Body
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Import model
fname = r'..\..\models\tbw.xbf'
bodies = Body.load_bodies(fname)

# Get bodies.
fuselage = bodies['Fuselage']
wing = bodies['Wing']
other_wing = bodies['Wing.2']
gear = bodies['Gear Pod']
other_gear = bodies['Gear Pod.1']
htail = bodies['Htail']
other_htail = bodies['Htail.3']
vtail = bodies['Vtail']
jury = bodies['Jury']
other_jury = bodies['Jury.4']
strut = bodies['Strut']
other_strut = bodies['Strut.5']

gui = Viewer()
for body in [fuselage, wing, other_wing, gear, other_gear, htail, other_htail,
             vtail, jury, other_jury, strut, other_strut]:
    body.set_color(0.5, 0.5, 0.5)
    body.set_transparency(0.5)
    gui.add(body)

# WING ------------------------------------------------------------------------
GroupAPI.create_group('wing group')

# Center wing structure will be based on the intersection between the wings
# and fuselage.
shape = IntersectShapes(fuselage.shape, wing.shape).shape

# Use the y-dimensions of the bounding box of the intersection curve.
bbox = BBox()
bbox.add_shape(shape)
ymax = bbox.ymax

# Center wing box
xz_plane = PlaneByAxes().plane
root_rib_pln = PlaneByAxes((0., ymax, 0.), 'xz').plane

# Front center spar
p1 = wing.sref.eval(0.15, 0.)
pln = PlaneByAxes(p1, 'yz').plane
fc_spar = SparBetweenShapes('fc spar', xz_plane, root_rib_pln, wing, pln).part

# Rear center spar
p1 = wing.sref.eval(0.70, 0.)
pln = PlaneByAxes(p1, 'yz').plane
rc_spar = SparBetweenShapes('rc spar', xz_plane, root_rib_pln, wing, pln).part

# Root rib
root_rib = RibByPoints('root rib', fc_spar.cref.p2, rc_spar.cref.p2, wing).part

# Strut rib where strut intersects.
p0 = strut.sref.eval(0.5, 0.)
p1 = strut.sref.eval(0.5, 1.)
strut_edge = EdgeByPoints(p0, p1).edge
bop = IntersectShapes(strut_edge, wing.shape)
v = bop.vertices[0]
p = v.point
strut_edge = EdgeByPoints(p0, p).edge
pln = PlaneByAxes(p, 'xz').plane
strut_rib = RibByShape('kink rib', pln, wing).part

# Tip rib
tip_rib = RibByParameters('tip rib', 0.15, 0.995, 0.70, 0.995, wing).part

# Inboard front spar
u = strut_rib.cref.local_to_global_u(0.15)
p2 = strut_rib.cref.eval(u)
pln = PlaneByIntersectingShapes(root_rib.shape, fc_spar.shape, p2).plane
inbd_fspar = SparByPoints('inbd fspar', root_rib.cref.p1, p2, wing, pln).part

# Inboard rear spar
u = strut_rib.cref.local_to_global_u(0.70)
p2 = strut_rib.cref.eval(u)
pln = PlaneByIntersectingShapes(root_rib.shape, rc_spar.shape, p2).plane
inbd_rspar = SparByPoints('inbd rspar', root_rib.cref.p2, p2, wing, pln).part

# Outboard front spar
pln = PlaneByIntersectingShapes(strut_rib.shape, inbd_fspar.shape,
                                tip_rib.cref.p1).plane
outbd_fspar = SparByPoints('outbd fspar', inbd_fspar.cref.p2, tip_rib.cref.p1,
                           wing, pln).part

# Outboard rear spar
pln = PlaneByIntersectingShapes(strut_rib.shape, inbd_rspar.shape,
                                tip_rib.cref.p2).plane
outbd_rspar = SparByPoints('outbd rspar', inbd_rspar.cref.p2, tip_rib.cref.p2,
                           wing, pln).part

# Jury rib where strut intersects
p0 = jury.sref.eval(0.5, 1.)
pln = PlaneByAxes(p0, 'xz').plane
jury_rib = RibBetweenShapes('jury rib', inbd_fspar.shape, inbd_rspar.shape,
                            wing, pln).part

# Inboard ribs
u2 = inbd_rspar.cref.invert(jury_rib.cref.p2)
inbd_ribs = RibsAlongCurveByDistance('inbd rib', inbd_rspar.cref, 70.,
                                     inbd_fspar.shape, inbd_rspar.shape, wing,
                                     u2=u2, d1=18, d2=-30).parts

# Middle ribs
mid_ribs = RibsAlongCurveByDistance('mid rib', inbd_rspar.cref, 70.,
                                    inbd_fspar.shape, inbd_rspar.shape, wing,
                                    u1=u2, d1=18, d2=-30).parts

# Outboard ribs
outbd_ribs = RibsAlongCurveByDistance('outbd rib', outbd_rspar.cref, 70.,
                                      outbd_fspar.shape, outbd_rspar.shape,
                                      wing, d1=18, d2=-30).parts

# Join and discard internal structure.
internal_parts = GroupAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

# Wing skin
wskin = SkinByBody('wing skin', wing).part

# Join the wing skin and internal structure
# wskin.fuse(*internal_parts)

# Discard faces touching reference surface.
# wskin.discard_by_dmin(wing.sref_shape, 1.)

# Fix skin since it's not a single shell anymore, but a compound of two
# shells (upper and lower skin).
# wskin.fix()
wskin.set_transparency(0.5)

# HTAIL -----------------------------------------------------------------------
# GroupAPI.create_group('htail group')
#
# fspar = SparByParameters('htail fspar', 0.15, 0.01, 0.15, 0.99, htail).part
# rspar = SparByParameters('htail rspar', 0.70, 0.01, 0.70, 0.99, htail).part
# root = RibByPoints('htail root rib', fspar.p1, rspar.p1, htail).part
# tip = RibByPoints('htail tip rib', fspar.p2, rspar.p2, htail).part
# ribs = RibsAlongCurveByDistance('htail rib', rspar.cref, 12, fspar.shape,
#                                 rspar.shape, htail, d1=12, d2=-12).parts
#
# internal_parts = GroupAPI.get_parts()
# FuseSurfacePartsByCref(internal_parts)
# DiscardByCref(internal_parts)
#
# skin = SkinByBody('htail skin', htail).part
# skin.fuse(*internal_parts)
# skin.discard_by_dmin(htail.sref_shape, 1.)
# skin.fix()
#
# gui.add(*GroupAPI.get_parts())

# VTAIL -----------------------------------------------------------------------
# GroupAPI.create_group('vtail group')
#
# u, v = ProjectPointToSurface(fspar.p1, vtail.sref).nearest_param
# mspar = SparByParameters('vtail mspar', u, 0.05, u, v, vtail).part
#
# u, v = ProjectPointToSurface(rspar.p1, vtail.sref).nearest_param
# rspar = SparByParameters('vtail rspar', u, 0.05, u, v, vtail).part
# fspar = SparByParameters('vtail fspar', 0.12, 0.05, 0.12, v, vtail).part
# RibByPoints('vtail root rib', fspar.p1, rspar.p1, vtail)
# RibByPoints('vtail tip rib', fspar.p2, rspar.p2, vtail)
# RibsAlongCurveByDistance('vtail rib', rspar.cref, 18., fspar.shape, rspar.shape,
#                          vtail, d1=18, d2=-72)
#
# internal_parts = GroupAPI.get_parts()
# FuseSurfacePartsByCref(internal_parts)
# DiscardByCref(internal_parts)
#
# skin = SkinByBody('vtail skin', vtail).part
# skin.fuse(*internal_parts)
# skin.discard_by_dmin(vtail.sref_shape, 1.)
# skin.fix()
#
# gui.add(*GroupAPI.get_parts())

# FUSELAGE --------------------------------------------------------------------
GroupAPI.create_group('fuselage group')

frame1 = FrameByPlane('frame 1', fc_spar.plane, fuselage, 4.).part
frame1.cut(wing)
frame1.cut(other_wing)

frame2 = FrameByPlane('frame 2', rc_spar.plane, fuselage, 4.).part
frame2.cut(wing)
frame2.cut(other_wing)

#
# pln = PlaneByAxes((60, 0, 0), 'yz').plane
# bh1 = BulkheadBySurface('bh 1', pln, fuselage).bulkhead
#
# pln = PlaneByAxes((180, 0, 0), 'yz').plane
# bh2 = BulkheadBySurface('bh 2', pln, fuselage).bulkhead
# p0 = fc_spar.p1
# pln = PlaneByAxes(p0, 'yz').plane
# bh3 = BulkheadBySurface('bh 3', pln, fuselage).bulkhead
#
# p0 = rc_spar.p1
# pln = PlaneByAxes(p0, 'yz').plane
# bh4 = BulkheadBySurface('bh 4', pln, fuselage).bulkhead
#
# p0 = fspar.p1
# pln = PlaneByAxes(p0, 'yz').plane
# bh5 = BulkheadBySurface('bh 5', pln, fuselage).bulkhead
#
# p0 = mspar.p1
# pln = PlaneByAxes(p0, 'yz').plane
# bh6 = BulkheadBySurface('bh 6', pln, fuselage).bulkhead
#
# p0 = rspar.p1
# pln = PlaneByAxes(p0, 'yz').plane
# bh7 = BulkheadBySurface('bh 7', pln, fuselage).bulkhead
#
# pln = PlaneByAxes((0, 0, -24), 'xy').plane
# floor = FloorBySurface('floor', pln, fuselage).floor
#
# shell = ShellByFaces(bh1.faces).shell
# fwd = HalfspaceByShape(shell, (-1e6, 0, 0)).solid
# floor.cut(fwd)
# shell = ShellByFaces(bh5.faces).shell
# aft = HalfspaceByShape(shell, (1e6, 0, 0)).solid
# floor.cut(aft)
#
# plns = [bh1.sref, bh2.sref, bh3.sref, bh4.sref, bh5.sref, bh6.sref, bh7.sref]
#
# indx = 1
# frames = []
# for pln1, pln2 in pairwise(plns):
#     builder = FramesBetweenPlanesByDistance('frame', pln1, pln2, 24., fuselage,
#                                             3.5, first_index=indx)
#     frames += builder.frames
#     indx = builder.next_index
#
# internal_parts = GroupAPI.get_parts()
#
# skin = SkinByBody('fuselage skin', fuselage).part
# skin.set_transparency(0.5)
#
# FuseSurfaceParts([skin], internal_parts)

# BEAMS -----------------------------------------------------------------------
GroupAPI.create_group('strut group')

beam1 = Beam1DByShape('beam 1', strut_edge).part
beam1.set_color(1, 0, 0)

p1 = jury.sref.eval(0.5, 0.5)
p2 = p1.copy()
v = VectorByArray((0, 0, 1.)).vector
ProjectPointToCurve(p1, beam1.cref, direction=v, update=True)
ProjectPointToCurve(p2, jury_rib.cref, direction=v, update=True)
beam2 = Beam1DByPoints('beam 2', p1, p2).part
beam2.set_color(1, 0, 0)
beam2.cut(wing)

# GEAR ------------------------------------------------------------------------
gear_group = GroupAPI.create_group('gear group')

# Spar at strut location
p = beam1.cref.eval(0.)
pln = PlaneByAxes(p, 'yz').plane
SparByShape('main gear spar', pln, gear)

x = SparByShape('gear spar 1', frame1.plane, gear).part
y = SparByShape('gear spar 2', frame2.plane, gear).part

frame1.cut(gear)
frame2.cut(gear)
frame1.cut(other_gear)
frame2.cut(other_gear)

# VIEW ------------------------------------------------------------------------
master = GroupAPI.get_master()
gui.add(master)
gui.start()
