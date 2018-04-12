from afem.adaptor import *
from afem.geometry import *
from afem.graphics import *
from afem.topology import *

gui = Viewer()

# Define some points to interpolate
pts = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]

# Interpolate the points
curve = NurbsCurveByInterp(pts).curve

# Create an edge from the curve
edge = EdgeByCurve(curve).edge
edge.set_color(1, 0, 0)

# Create an adaptor curve from the edge
adp_crv = EdgeAdaptorCurve.by_edge(edge)

# Adaptors work in some geometry tools like projection
p1 = Point(8, 8, 0)
proj = ProjectPointToCurve(p1, adp_crv)

# Create an edge to shown projection
e1 = EdgeByPoints(p1, proj.nearest_point).edge
e1.set_color(0, 0, 1)

# View the results
gui.add(edge, p1, e1)
gui.view_top()
gui.start()

# Create a face by dragging the edge
face = FaceByDrag(edge, (0, 0, 10)).face
face.set_color(0, 1, 0)

# Create an adaptor surface from the face
adp_srf = FaceAdaptorSurface.by_face(face)

# Project a point to the adaptor surface
p2 = Point(8, 8, 5)
proj = ProjectPointToSurface(p2, adp_srf)

# Create an edge to shown projection
e2 = EdgeByPoints(p2, proj.nearest_point).edge
e2.set_color(0, 0, 1)

# View the results
gui.add(face, p2, e2)
gui.view_iso()
gui.start()

# Adaptors have similar interfaces to curves and surfaces
p3 = adp_crv.eval(10)
p4 = adp_srf.eval(10, -5)

# View the points
gui.add(p3, p4)
gui.start()

# Create points along an adaptor curve
tool = PointsAlongCurveByNumber(adp_crv, 10)

# View the points
gui.clear()
gui.add(edge, *tool.points)
gui.view_top()
gui.start()

# Define some new points to interpolate
pts = [(10, 0, 0), (15, 5, 0), (20, 0, 0)]

# Interpolate the points
curve2 = NurbsCurveByInterp(pts).curve

# Create an edge from the curve
edge2 = EdgeByCurve(curve2).edge
edge.set_color(1, 0, 0)

# Build a wire from the edges
wire = WiresByConnectedEdges([edge, edge2]).wires[0]
wire.set_color(1, 0, 0)

# Create an adaptor curve from a wire. The static method "to_adaptor" can be
# used for convenience to convert entities to an adaptor if possible.
adp_crv2 = AdaptorCurve.to_adaptor(wire)

# Create points along the wire adaptor curve
tool = PointsAlongCurveByNumber(adp_crv2, 20)

# View the points
gui.clear()
gui.add(wire, *tool.points)
gui.view_top()
gui.start()
