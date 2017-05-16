from __future__ import print_function

from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import AssemblyData, CreatePart, PartTools

# Inputs
fname = r'..\models\N2A_nosplit.stp'

# Import model
ImportVSP.step_file(fname)
wing = ImportVSP.get_body('Wing_Body')
other_wing = ImportVSP.get_body('Wing_Body.1')
vtail = ImportVSP.get_body('Vertical_Tails')
other_vtail = ImportVSP.get_body('Vertical_Tails.2')
for name in ImportVSP.get_bodies():
    body = ImportVSP.get_body(name)
    body.set_transparency(0.5)
    body.set_color(0.5, 0.5, 0.5)
    Viewer.add_items(body)

# Outboard wing structure

# Root rib
outbd_root_rib = CreatePart.rib.by_parameters('outbd root rib', wing,
                                              0.15, 0.5, 0.65, 0.5)

# Tip rib
tip_rib = CreatePart.rib.by_parameters('tip rib', wing,
                                       0.15, 0.995, 0.65, 0.995)

# Outbd front spar
outbd_fspar = CreatePart.spar.by_points('outbd front spar', wing,
                                        outbd_root_rib.p1, tip_rib.p1)

# Outbd rear spar
outbd_rspar = CreatePart.spar.by_points('outbd rear spar', wing,
                                        outbd_root_rib.p2, tip_rib.p2)

# Outbd ribs
outbd_ribs = CreatePart.rib.along_curve('outbd rib', wing, outbd_rspar.cref,
                                        outbd_fspar.sref, outbd_rspar.sref,
                                        30., s1=30., s2=-30.)

internal_parts = AssemblyData.get_parts()

PartTools.fuse_wing_parts(internal_parts)
PartTools.discard_faces(internal_parts)

Viewer.add_items(*internal_parts)

Viewer.show()
