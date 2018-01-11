from numpy import linspace

from afem.config import Settings
from afem.exchange import StepImport, ImportVSP
from afem.graphics import Viewer
from afem.misc.check import pairwise
from afem.topology import *

Settings.log_to_console()

# Set units to inch
Settings.set_units('in')

# OpenVSP 3.5.0 ----------------------------------------------------------------

fn = r'../models/777-200LR_nosplit.stp'
step = StepImport()
step.read(fn)

ImportVSP.step_file(fn, False)
wing = ImportVSP.get_body('Wing')
vsp_srf_35 = wing.get_metadata('vsp surface')

v = Viewer()
v.view.set_white_background()
v.view.display_shape(step.shape, (0.5, 0.5, 0.5))
v.start()
v.clear()

c0_shape = DivideC0Shape(step.shape).shape

for f in ExploreShape.get_faces(c0_shape):
    v.view.display_shape(f)
v.start()
v.clear()

# OpenVSP 3.15.0 ---------------------------------------------------------------

fn = r'../models/777-200LR_nosplit_vsp315.stp'
step = StepImport()
step.read(fn)

ImportVSP.step_file(fn, False)
wing = ImportVSP.get_body('Wing')
vsp_srf_315 = wing.get_metadata('vsp surface')

v.view.display_shape(step.shape, (0.5, 0.5, 0.5))
v.start()
v.clear()

c0_shape = DivideC0Shape(step.shape).shape

for f in ExploreShape.get_faces(c0_shape):
    v.view.display_shape(f)
v.start()
v.view.remove_all()

# First derivatives ------------------------------------------------------------
u_35 = vsp_srf_35.local_to_global_param('u', 0.5)
u_315 = vsp_srf_315.local_to_global_param('u', 0.5)

# Make sure curves are at same location
u_iso_35 = vsp_srf_35.u_iso(u_35)
u_iso_315 = vsp_srf_315.u_iso(u_315)
v.add(u_iso_35, u_iso_315, vsp_srf_35)
v.start()

# Evaluate v-direction first derivative at 100 points between each knot value.
# Record the v-parameter and magnitude to a csv file.
n = 100
fout = open('dv_35.csv', 'w')
for v1, v2 in pairwise(vsp_srf_35.vknots):
    for v in linspace(v1, v2, n):
        dv = vsp_srf_35.deriv(u_35, v, 0, 1)
        fout.write('{}, {}\n'.format(v, dv.mag))

fout = open('dv_315.csv', 'w')
for v1, v2 in pairwise(vsp_srf_315.vknots):
    for v in linspace(v1, v2, n):
        dv = vsp_srf_315.deriv(u_315, v, 0, 1)
        fout.write('{}, {}\n'.format(v, dv.mag))
