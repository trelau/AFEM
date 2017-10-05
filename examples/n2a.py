from afem.config import Settings
from afem.geometry import *
from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.structure import *
from afem.topology import *

Settings.log_to_console(True)

# Inputs
fname = r'..\models\N2A_nosplit.stp'
# Length of flight deck from nose.
fd_length = 10. * 12.
# Cabin length as absolute value (> 1) or percent root chord ( < 1)
cabin_length = 0.7
# Cabin width
cabin_width = 50. * 12.
# Bay width
bay_width = 8. * 12.

# Import model
ImportVSP.step_file(fname)
wing = ImportVSP.get_body('Wing_Body')
other_wing = ImportVSP.get_body('Wing_Body.1')
vtail = ImportVSP.get_body('Vertical_Tails')
other_vtail = ImportVSP.get_body('Vertical_Tails.2')
for body in [wing, other_wing, vtail, other_vtail]:
    body.set_transparency(0.5)
    body.set_color(0.5, 0.5, 0.5)

# Construction geometry
root_chord = wing.isocurve(v=0.)

AssemblyAPI.create_assy('wing assy')

# Centerbody structure
pln = PlaneByAxes((fd_length, 0, 0), 'yz').plane
fd_bh = SparBySurface('flight deck bh', pln, wing).spar

# Rear cabin bulkhead
if cabin_length < 1:
    ds = cabin_length * root_chord.length
    p0 = PointFromParameter(root_chord, root_chord.u1, ds).point
    pln = PlaneByAxes(p0, 'yz').plane
else:
    p0 = Point(fd_length, cabin_length, 0, 0)
    pln = PlaneByAxes(p0, 'yz').plane
rear_cabin_bh = SparBySurface('rear cabin bulkhead', pln, wing).spar

# Cabin side-of-body
p0 = Point(p0.x, cabin_width / 2., 0)
pln = PlaneByAxes(p0, 'xz').plane
cb_sob = RibBySurface('centerbody sob', pln, wing).rib

# Trim cref of rear bulkhead
rear_cabin_bh.trim_u2(cb_sob)

# Cabin bay wall
p0 = Point(p0.x, cabin_width / 4., 0)
pln = PlaneByAxes(p0, 'xz').plane
cb_wall = RibBySurface('centerbody wall', pln, wing).rib

# Front spar segments
p1 = fd_bh.point_from_parameter(0.75, is_rel=True)
p2 = cb_wall.point_from_parameter(0.05, is_rel=True)
cb_fspar1 = SparByPoints('cb spar 1', p1, p2, wing).spar

p2 = cb_sob.point_from_parameter(0.05, is_rel=True)
pln = PlaneByIntersectingShapes(cb_wall, cb_fspar1, p2).plane
cb_fspar2 = SparByPoints('cb fspar 2', cb_fspar1.p2, p2, wing, pln).spar

# Outboard wing structure

# Root rib
root_rib = RibByParameters('ob root rib', 0.15, 0.5, 0.7, 0.5, wing).rib

# Tip rib
tip_rib = RibByParameters('ob tip rib', 0.15, 0.995, 0.7, 0.995, wing).rib

# Outbd front spar
ob_fspar = SparByPoints('ob fspar', root_rib.p1, tip_rib.p1, wing).spar

# Outbd rear spar
ob_rspar = SparByPoints('ob rspar', root_rib.p2, tip_rib.p2, wing).spar

# Outbd ribs
RibsAlongCurveByDistance('ob rib', ob_rspar.cref, 30, ob_fspar,
                         ob_rspar, wing, d1=30, d2=-30)

# Trap wing (need to redo and align spars given intersection of parts)
fspar = SparByPoints('trap fspar', cb_fspar2.p2, ob_fspar.p1, wing).spar
rspar = SparByPoints('trap rspar', rear_cabin_bh.p2, ob_rspar.p1, wing).spar

RibsBetweenPlanesByDistance('trap rib', cb_sob.plane, root_rib.plane, 30,
                            fspar, rspar, wing, 30, -30)

# Rear centerbody
p0 = root_chord.eval(0.85)
pln = PlaneByAxes(p0, 'yz').plane
SparBySurface('rc spar', pln, wing)

# Spars and cb intersections
p1 = cb_fspar1.p2
pln = PlaneByAxes(p1, 'yz').plane
spar1 = SparBySurface('cb spar 1', pln, wing).spar
spar1.trim_u2(cb_wall)

p1 = cb_fspar2.p2
ProjectPointToCurve(p1, root_chord, update=True)
pln = PlaneByAxes(p1, 'yz').plane
spar2 = SparBySurface('cb spar 2', pln, wing).spar
spar2.trim_u2(cb_sob)

p0 = rear_cabin_bh.plane.eval()
p0.x -= 144
pln = PlaneByAxes(p0, 'yz').plane
spar3 = SparBySurface('cb spar 3', pln, wing).spar
spar3.trim_u2(cb_sob)

# Cut with an imaginary floor for now
pln = PlaneByAxes((0, 0, 0), axes='xy').plane
solid = HalfspaceBySurface(pln, (0, 0, 100)).solid
spar1.cut(solid)
spar2.cut(solid)
spar3.cut(solid)

# Fuse
internal_parts = AssemblyAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

skin = SkinByBody('skin', wing).skin
skin.fuse(*internal_parts)
skin.set_transparency(0.5)

# Vtail structure
AssemblyAPI.create_assy('vtail assy')

fspar = SparByParameters('vtail fspar', 0.15, 0.01, 0.15, 0.99, vtail).spar
rspar = SparByParameters('vtail rspar', 0.70, 0.01, 0.70, 0.99, vtail).spar
RibByPoints('vtail root rib', fspar.p1, rspar.p1, vtail)
RibByPoints('vtail tip rib', fspar.p2, rspar.p2, vtail)
RibsAlongCurveByDistance('vtail rib', rspar.cref, 18, fspar, rspar,
                         vtail, d1=18, d2=-30)

internal_parts = AssemblyAPI.get_parts()
FuseSurfacePartsByCref(internal_parts)
DiscardByCref(internal_parts)

skin = SkinByBody('vtail skin', vtail).skin
skin.fuse(*internal_parts)
skin.discard_by_dmin(vtail.sref_shape, 1.)
skin.fix()

# Mirror and plot
parts = AssemblyAPI.get_master().get_parts()
xz_pln = PlaneByAxes().plane
for part in parts:
    part.set_mirror(xz_pln)

Viewer.add(*AssemblyAPI.get_master().get_parts())
Viewer.show()
