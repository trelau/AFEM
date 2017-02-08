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

fuselage.set_transparency(0.5)
wing.set_transparency(0.5)

# Fwd bulkhead at 25% wing chord
p0 = wing.eval(0.25, 0.)
pln1 = CreateGeom.plane_by_axes(p0, 'yz')
fwd_cut = ShapeTools.box_from_plane(pln1, 500, 500, -500)

# Rear bulkhead at 65% wing chord
p0 = wing.eval(0.65, 0.)
pln2 = CreateGeom.plane_by_axes(p0, 'yz')
aft_cut = ShapeTools.box_from_plane(pln2, 500, 500, 500)

# Cut for lower wing
pref = wing.eval(0.5, 0.5)
lower_cut = ShapeTools.make_halfspace(wing.sref, [0, 0, -1000])

# Find intersection of fuselage and wing
shape = fuselage.section(wing)

# Cut away parts of intersection curve to get only upper wing intersection
shape = ShapeTools.bcut(shape, fwd_cut)
shape = ShapeTools.bcut(shape, aft_cut)
shape = ShapeTools.bcut(shape, lower_cut)

# Make face in z-direction
vz = CreateGeom.vector([0, 0, -100])
shape = ShapeTools.make_prism(shape, vz)

# Build the rib
rib = CreatePart.surface_part('rib', shape, wing)
rib.build()

Viewer.add_items(fuselage, wing, rib)
Viewer.show()
