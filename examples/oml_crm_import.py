from afem.exchange import brep
from afem.geometry import NurbsCurveByPoints, NurbsSurfaceByInterp
from afem.graphics import Viewer
from afem.oml import Body

fn1 = r'..\models\uCRM\fuselage.brep'
fn2 = r'..\models\uCRM\lhs_wing.brep'
fn3 = r'..\models\uCRM\rhs_wing.brep'

shape1 = brep.read_brep(fn1)
shape2 = brep.read_brep(fn2)
shape3 = brep.read_brep(fn3)

fuselage = Body(shape1, 'fuselage')
lhs_wing = Body(shape2, 'lhs wing')
rhs_wing = Body(shape3, 'rhs wing')

fuselage.set_transparency(0.5)
lhs_wing.set_transparency(0.5)
rhs_wing.set_transparency(0.5)

# Manually build wing reference surfaces
p1 = [23.06578, 0., 4.423368]
p2 = [36.52671, 0., 2.857832]
c1 = NurbsCurveByPoints([p1, p2]).curve

p1 = [25.217943, 3.054401, 4.46611]
p2 = [37.050119, 3.054681, 3.562711]
c2 = NurbsCurveByPoints([p1, p2]).curve

p1 = [31.126975, 10.877183, 4.513347]
p2 = [38.382418, 10.890519, 4.304047]
c3 = NurbsCurveByPoints([p1, p2]).curve

p1 = [45.164217, 29.445051, 4.63929]
p2 = [47.890914, 29.45361, 4.622345]
c4 = NurbsCurveByPoints([p1, p2]).curve

sref = NurbsSurfaceByInterp([c1, c2, c3, c4], 1).surface
rhs_wing.set_sref(sref)

p1 = [23.06578, 0., 4.423368]
p2 = [36.52671, 0., 2.857832]
c1 = NurbsCurveByPoints([p1, p2]).curve

p1 = [25.217943, -3.054401, 4.46611]
p2 = [37.050119, -3.054681, 3.562711]
c2 = NurbsCurveByPoints([p1, p2]).curve

p1 = [31.126975, -10.877183, 4.513347]
p2 = [38.382418, -10.890519, 4.304047]
c3 = NurbsCurveByPoints([p1, p2]).curve

p1 = [45.164217, -29.445051, 4.63929]
p2 = [47.890914, -29.45361, 4.622345]
c4 = NurbsCurveByPoints([p1, p2]).curve

sref = NurbsSurfaceByInterp([c1, c2, c3, c4], 1).surface
lhs_wing.set_sref(sref)

gui = Viewer()
for shape in [fuselage, lhs_wing, rhs_wing]:
    gui.add(shape)
gui.start()
