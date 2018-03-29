from afem.geometry import *
from afem.graphics import Viewer

gui = Viewer()

# Create a point directly from the entity. Default is (0, 0, 0).
p1 = Point()

# Create a point by array-like
p2 = PointByArray([5, 0, 5]).point

# Create a point by x-, y-, and z-coordinates.
p3 = PointByXYZ(10, 0, 0).point

# Interpolate the points with a curve
c1 = NurbsCurveByInterp([p1, p2, p3]).curve

gui.add(p1, p2, p3, c1)
gui.start()

# Copy curve and translate
c2 = c1.copy()
c2.translate((0, 10, 0))

gui.add(c2)
gui.start()

# Copy and translate again
c3 = c2.copy()
c3.translate((0, 10, 10))

gui.add(c3)
gui.start()

# Approximate a surface
s1 = NurbsSurfaceByApprox([c1, c2, c3]).surface

gui.add(s1)
gui.start()

# Extract an iso-curve
c4 = s1.u_iso(10.)

gui.add(c4)
gui.start()

# Create points along the curve
pnts = PointsAlongCurveByDistance(c4, 1.).points

gui.add(*pnts)
gui.start()

# Extract iso-curve
c5 = s1.v_iso(0.5)

gui.add(c5)
gui.start()

# Intersect two curves
cci = IntersectCurveCurve(c4, c5)

gui.clear()
gui.add(c4, c5, s1, *cci.points)
gui.start()

# Define a plane along a curve
pln = PlaneFromParameter(c4, 0., 2.).plane

# Intersect a surface and a plane
ssi = IntersectSurfaceSurface(s1, pln)

gui.add(s1, *ssi.curves)
gui.start()

# Project a point to a surface
p4 = pln.eval(5, 5)
proj = ProjectPointToSurface(p4, s1)
line = NurbsCurveByInterp([p4, proj.nearest_point]).curve

gui.add(p4, proj.nearest_point, line)
gui.start()
