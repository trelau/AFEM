from afem.config import Settings
from afem.exchange import *
from afem.graphics import *
from afem.topology import *

Settings.log_to_console()

v = Viewer()

# # Default import
# fname = r'..\models\777-200LR.stp'
# fname = r'..\models\777-200LR_nosplit.stp'
# ImportVSP.step_file(fname)
# wing = ImportVSP.get_body('Wing')
# v.add(wing)
# v.start()
# v.clear()
#
# split_shape = DivideC0Shape(wing.solid).shape
# v.add(split_shape)
# v.start()
# v.clear()

# Import at try BSpline restriction (auto refit with BSplines). Note that there
# is some surface error but the chord-wise C0 splits are gone. This refits the
# surfaces after initial solid construction. Could take a while to run,
# especially for non-split wings and newer versions of OpenVSP.
# fname = r'..\models\777-200LR.stp'
# fname = r'..\models\777-200LR_nosplit.stp'
# ImportVSP.step_file(fname, bspline_restrict=True, tol=0.01)
# wing = ImportVSP.get_body('Wing')
# v.add(wing)
# v.start()
# v.clear()
#
# split_shape = DivideC0Shape(wing.solid).shape
# v.add(split_shape)
# v.start()
# v.clear()

# Try to reloft a wing using curves at each spanwise section. Very dependent on
# tessellation tolerance and doesn't support blended wing sections and forces
# flat end caps. The lower the tolerance the longer it takes to run, but too
# high a tolerance causes self-intersecting results if the TE is thin. Just
# doesn't seem to be a general/robust solution.
fname = r'..\models\777-200LR_nosplit.stp'
ImportVSP.step_file(fname, reloft=True, tol=0.01)
# tol=0.005 causes issues...
# ImportVSP.step_file(fname, reloft=True, tol=0.005)
wing = ImportVSP.get_body('Wing')
v.add(wing)
v.start()
v.clear()

split_shape = DivideC0Shape(wing.solid).shape
v.add(split_shape)
v.start()
v.clear()
