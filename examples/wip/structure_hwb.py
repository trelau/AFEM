from afem.config import Settings
from afem.exchange import ImportVSP
from afem.geometry import *
from afem.graphics import Viewer
from afem.structure import *
from afem.topology import *

Settings.log_to_console()

# Inputs
fname = r'..\..\models\hwb_nosplit.stp'
# Length of flight deck from nose.
fd_length = 10. * 12.
# Cabin length as absolute value (> 1) or percent root chord ( < 1)
cabin_length = 0.7
# Cabin width
cabin_width = 50. * 12.
# Bay width
bay_width = 8. * 12.

# Import model
vsp_import = ImportVSP()
vsp_import.import_step(fname)
wing = vsp_import['Wing_Body']
other_wing = vsp_import['Wing_Body.1']
vtail = vsp_import['Vertical_Tails']
other_vtail = vsp_import['Vertical_Tails.2']
for body in [wing, other_wing, vtail, other_vtail]:
    body.set_transparency(0.5)
    body.set_color(0.5, 0.5, 0.5)

# Construction geometry
root_chord = wing.sref.v_iso(0.)

GroupAPI.create_group('wing group')

# Centerbody structure
pln = PlaneByAxes((fd_length, 0, 0), 'yz').plane
fd_bh = SparByShape('flight deck bh', pln, wing).part

# Rear cabin bulkhead
if cabin_length < 1:
    ds = cabin_length * root_chord.length
    p0 = PointFromParameter(root_chord, root_chord.u1, ds).point
    pln = PlaneByAxes(p0, 'yz').plane
else:
    p0 = Point(fd_length, cabin_length, 0, 0)
    pln = PlaneByAxes(p0, 'yz').plane
# TODO The script fails here due to Boolean operation...
rear_cabin_bh = SparByShape('rear cabin bulkhead', pln, wing).part

# Cabin side-of-body
p0 = Point(p0.x, cabin_width / 2., 0)
pln = PlaneByAxes(p0, 'xz').plane
cb_sob = RibByShape('centerbody sob', pln, wing).part

# Trim cref of rear bulkhead
rear_cabin_bh.trim_u2(cb_sob.shape)

# Cabin bay wall
p0 = Point(p0.x, cabin_width / 4., 0)
pln = PlaneByAxes(p0, 'xz').plane
cb_wall = RibByShape('centerbody wall', pln, wing).part

# Front spar segments
p1 = fd_bh.point_from_parameter(0.75, is_rel=True)
p2 = cb_wall.point_from_parameter(0.05, is_rel=True)
cb_fspar1 = SparByPoints('cb spar 1', p1, p2, wing).part

p2 = cb_sob.point_from_parameter(0.05, is_rel=True)
pln = PlaneByIntersectingShapes(cb_wall.shape, cb_fspar1.shape, p2).plane
cb_fspar2 = SparByPoints('cb fspar 2', cb_fspar1.cref.p2, p2, wing, pln).part

# Outboard wing structure

# Root rib
root_rib = RibByParameters('ob root rib', 0.15, 0.5, 0.7, 0.5, wing).part

# Tip rib
tip_rib = RibByParameters('ob tip rib', 0.15, 0.995, 0.7, 0.995, wing).part

# Outbd front spar
ob_fspar = SparByPoints('ob fspar', root_rib.cref.p1, tip_rib.cref.p1,
                        wing).part

# Outbd rear spar
ob_rspar = SparByPoints('ob rspar', root_rib.cref.p2, tip_rib.cref.p2,
                        wing).part

# Outbd ribs
RibsAlongCurveByDistance('ob rib', ob_rspar.cref, 30, ob_fspar.shape,
                         ob_rspar.shape, wing, d1=30, d2=-30)

# Trap wing (need to redo and align spars given intersection of parts)
fspar = SparByPoints('trap fspar', cb_fspar2.cref.p2, ob_fspar.cref.p1,
                     wing).part
rspar = SparByPoints('trap rspar', rear_cabin_bh.cref.p2, ob_rspar.cref.p1,
                     wing).part

RibsBetweenPlanesByDistance('trap rib', cb_sob.plane, root_rib.plane, 30,
                            fspar.shape, rspar.shape, wing, 30, -30)

# Rear centerbody
p0 = root_chord.eval(0.85)
pln = PlaneByAxes(p0, 'yz').plane
SparByShape('rc spar', pln, wing)

# Spars and cb intersections
p1 = cb_fspar1.cref.p2
pln = PlaneByAxes(p1, 'yz').plane
spar1 = SparByShape('cb spar 1', pln, wing).part
spar1.trim_u2(cb_wall.shape)

p1 = cb_fspar2.cref.p2
ProjectPointToCurve(p1, root_chord, update=True)
pln = PlaneByAxes(p1, 'yz').plane
spar2 = SparByShape('cb spar 2', pln, wing).part
spar2.trim_u2(cb_sob.shape)

p0 = rear_cabin_bh.plane.eval()
p0.x -= 144
pln = PlaneByAxes(p0, 'yz').plane
spar3 = SparByShape('cb spar 3', pln, wing).part
spar3.trim_u2(cb_sob.shape)

# Cut with an imaginary floor for now
pln = PlaneByAxes((0, 0, 0), axes='xy').plane
solid = HalfspaceBySurface(pln, (0, 0, 100)).solid
spar1.cut(solid)
spar2.cut(solid)
spar3.cut(solid)

# Fuse
internal_parts = GroupAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

skin = SkinByBody('skin', wing).part
skin.fuse(*internal_parts)
skin.set_transparency(0.5)

# Vtail structure
GroupAPI.create_group('vtail group')

fspar = SparByParameters('vtail fspar', 0.15, 0.01, 0.15, 0.99, vtail).part
rspar = SparByParameters('vtail rspar', 0.70, 0.01, 0.70, 0.99, vtail).part
RibByPoints('vtail root rib', fspar.cref.p1, rspar.cref.p1, vtail)
RibByPoints('vtail tip rib', fspar.cref.p2, rspar.cref.p2, vtail)
RibsAlongCurveByDistance('vtail rib', rspar.cref, 18, fspar.shape, rspar.shape,
                         vtail, d1=18, d2=-30)

internal_parts = GroupAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

skin = SkinByBody('vtail skin', vtail).part
skin.fuse(*internal_parts)
skin.discard_by_dmin(vtail.sref_shape, 1.)
skin.fix()

# View
skin.set_transparency(0.5)
gui = Viewer()
gui.add(GroupAPI.get_master())
gui.start()
