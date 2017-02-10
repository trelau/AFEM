from __future__ import print_function

import time

from asap.graphics import Viewer
from asap.io import ImportVSP
from asap.structure import CreatePart

# Import model
fn = r'.\test_io\777-200LR_mod_vsp350_sref.stp'
ImportVSP.step_file(fn)

# Get components.
wing = ImportVSP.get_body('Wing')
fuselage = ImportVSP.get_body('Fuselage')

start = time.time()
# Spars
fspar = CreatePart.spar.by_parameters('fspar', wing, 0.15, 0.05, 0.15, 0.925)
rspar = CreatePart.spar.by_parameters('rspar', wing, 0.65, 0.05, 0.65, 0.925)

fspar.cut_hole(12. * 20, 10.)
fspar.cut_hole(12. * 15, 15.)

# Bulkheads
bh1 = CreatePart.bulkhead.by_sref('bh1', fuselage, fspar.rshape)
bh2 = CreatePart.bulkhead.by_sref('bh2', fuselage, rspar.rshape)

# Join
fspar.fuse(bh1)
rspar.fuse(bh2)

for part in [fspar, rspar, bh1, bh2]:
    part.mesh(4., quad_dominated=False)
Viewer.add_meshes(bh1, bh2, fspar, rspar)
Viewer.show_mesh()

Viewer.add_items(fspar, rspar, bh1, bh2)
Viewer.show()

# Now test with sewing
fspar = CreatePart.spar.by_parameters('fspar', wing, 0.15, 0.05, 0.15, 0.925)
rspar = CreatePart.spar.by_parameters('rspar', wing, 0.65, 0.05, 0.65, 0.925)

bh1 = CreatePart.bulkhead.by_sref('bh1', fuselage, fspar.rshape)
bh2 = CreatePart.bulkhead.by_sref('bh2', fuselage, rspar.rshape)

# Cut bulkheads with spars
bh1.cut(fspar)
bh2.cut(rspar)

Viewer.add_items(bh1, bh2, fspar, rspar)
Viewer.show()

# Shared edges should be empty
print(bh1.shared_edges(fspar))
print(bh2.shared_edges(rspar))

# Now sew
fspar.sew(bh1)
rspar.sew(bh2)

Viewer.add_items(bh1, bh2, fspar, rspar)
Viewer.show()

for part in [fspar, rspar]:
    part.mesh(4., quad_dominated=False)
for part in [bh1, bh2]:
    part.mesh(8., quad_dominated=False)
Viewer.add_meshes(bh1, bh2, fspar, rspar)
Viewer.show_mesh()
