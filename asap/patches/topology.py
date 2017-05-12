from OCC.BRep import BRep_Tool_IsClosed, BRep_Tool_Parameter, BRep_Tool_Pnt
from OCC.BRepAdaptor import BRepAdaptor_CompCurve, BRepAdaptor_Curve, \
    BRepAdaptor_Surface
from OCC.BRepCheck import BRepCheck_Analyzer
from OCC.BRepClass3d import brepclass3d_OuterShell
from OCC.BRepPrimAPI import BRepPrimAPI_MakeHalfSpace
from OCC.BRepTools import breptools_OuterWire
from OCC.ShapeAnalysis import ShapeAnalysis_Edge
from OCC.TopAbs import TopAbs_COMPOUND, TopAbs_COMPSOLID, TopAbs_EDGE, \
    TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE
from OCC.TopoDS import TopoDS_Edge, TopoDS_Face, TopoDS_Shape, TopoDS_Shell, \
    TopoDS_Solid, TopoDS_Vertex, TopoDS_Wire

from ..geometry.points import Point
from ..topology import ShapeTools

_analysis_edge = ShapeAnalysis_Edge()


def _to_shape(self):
    """
    Convert to shape.
    """
    return ShapeTools.to_shape(self)


def _is_vertex(self):
    return self.ShapeType() == TopAbs_VERTEX


def _is_edge(self):
    return self.ShapeType() == TopAbs_EDGE


def _is_wire(self):
    return self.ShapeType() == TopAbs_WIRE


def _is_face(self):
    return self.ShapeType() == TopAbs_FACE


def _is_shell(self):
    return self.ShapeType() == TopAbs_SHELL


def _is_solid(self):
    return self.ShapeType() == TopAbs_SOLID


def _is_compsolid(self):
    return self.ShapeType() == TopAbs_COMPSOLID


def _is_compound(self):
    return self.ShapeType() == TopAbs_COMPOUND


def _get_vertices(self):
    """
    Get vertices of shape.
    """
    return ShapeTools.get_vertices(self)


def _get_edges(self):
    """
    Get edges of shape.
    """
    return ShapeTools.get_edges(self)


def _get_wires(self):
    """
    Get wires of shape.
    """
    return ShapeTools.get_wires(self)


def _get_faces(self):
    """
    Get faces of shape.
    """
    return ShapeTools.get_faces(self)


def _get_shells(self):
    """
    Get shell of shape.
    """
    return ShapeTools.get_shells(self)


def _get_solids(self):
    """
    Get solids of shape.
    """
    return ShapeTools.get_solids(self)


def _get_compounds(self):
    """
    Get compounds of shape.
    """
    return ShapeTools.get_compounds(self)


def _is_closed(self):
    """
    Is shape closed.
    """
    return BRep_Tool_IsClosed(self)


def _is_valid(self):
    """
    Check validity of shape.
    """
    check = BRepCheck_Analyzer(self, True)
    return check.IsValid()


def _vertex_pnt(self):
    """
    Point at vertex.
    """
    gp_pnt = BRep_Tool_Pnt(self)
    return Point(gp_pnt.XYZ())


def _edge_eval(self, u):
    """
    Evaluate point on edge.
    """
    adp_crv = BRepAdaptor_Curve(self)
    p = Point()
    adp_crv.D0(u, p)
    return p


def _wire_eval(self, u):
    """
    Evaluate point on wire.
    """
    adp_crv = BRepAdaptor_CompCurve(self)
    p = Point()
    adp_crv.D0(u, p)
    return p


def _face_eval(self, u, v):
    """
    Evaluate point on face.
    """
    adp_srf = BRepAdaptor_Surface(self)
    p = Point()
    adp_srf.D0(u, v, p)
    return p


def _edge_v1(self):
    """
    First vertex of edge.
    """
    return _analysis_edge.FirstVertex(self)


def _edge_v2(self):
    """
    Last vertex of edge.
    """
    return _analysis_edge.LastVertex(self)


def _edge_u1(self):
    """
    First parameter of edge.
    """
    v1 = _edge_v1(self)
    return BRep_Tool_Parameter(v1, self)


def _edge_u2(self):
    """
    Last parameter of edge.
    """
    v2 = _edge_v2(self)
    return BRep_Tool_Parameter(v2, self)


def _edge_p1(self):
    """
    First point of edge. 
    """
    v1 = _edge_v1(self)
    return _vertex_pnt(v1)


def _edge_p2(self):
    """
    Last point of edge. 
    """
    v2 = _edge_v2(self)
    return _vertex_pnt(v2)


def _ordered_edges(self):
    """
    List of wire edges from wire explorer.
    """


def _outer_wire(self):
    return breptools_OuterWire(self)


def _outer_shell(self):
    return brepclass3d_OuterShell(self)


def _fuse(self, other, rtype=''):
    return ShapeTools.bfuse(self, other, rtype)


def _common(self, other, rtype=''):
    return ShapeTools.bcommon(self, other, rtype)


def _section(self, other, rtype=''):
    return ShapeTools.bsection(self, other, rtype)


def _cut(self, other, rtype=''):
    return ShapeTools.bcut(self, other, rtype)


def _make_halfspace(self, pref):
    """
    Make a halfspace of the face or shell.
    """
    return BRepPrimAPI_MakeHalfSpace(self, pref).Solid()


def _mass(self):
    return ShapeTools.shape_volume(self)


def _distance(self, other):
    return ShapeTools.min_distance(self, other)[0]


TopoDS_Shape.shape = property(_to_shape)
TopoDS_Shape.is_closed = property(_is_closed)
TopoDS_Shape.is_valid = property(_is_valid)
TopoDS_Shape.is_vertex = property(_is_vertex)
TopoDS_Shape.is_edge = property(_is_edge)
TopoDS_Shape.is_wire = property(_is_wire)
TopoDS_Shape.is_face = property(_is_face)
TopoDS_Shape.is_shell = property(_is_shell)
TopoDS_Shape.is_solid = property(_is_solid)
TopoDS_Shape.is_compsolid = property(_is_compsolid)
TopoDS_Shape.is_compound = property(_is_compound)

TopoDS_Shape.vertices = property(_get_vertices)
TopoDS_Shape.edges = property(_get_edges)
TopoDS_Shape.wires = property(_get_wires)
TopoDS_Shape.faces = property(_get_faces)
TopoDS_Shape.shells = property(_get_shells)
TopoDS_Shape.solids = property(_get_solids)
TopoDS_Shape.compounds = property(_get_compounds)

TopoDS_Shape.fuse = _fuse
TopoDS_Shape.common = _common
TopoDS_Shape.section = _section
TopoDS_Shape.cut = _cut
TopoDS_Shape.distance = _distance

TopoDS_Vertex.pnt = property(_vertex_pnt)

TopoDS_Edge.eval = _edge_eval
TopoDS_Edge.v1 = property(_edge_v1)
TopoDS_Edge.v2 = property(_edge_v2)
TopoDS_Edge.u1 = property(_edge_u1)
TopoDS_Edge.u2 = property(_edge_u2)
TopoDS_Edge.p1 = property(_edge_p1)
TopoDS_Edge.p2 = property(_edge_p2)
TopoDS_Edge.length = property(_mass)

TopoDS_Wire.eval = _wire_eval
TopoDS_Wire.length = property(_mass)

TopoDS_Face.eval = _face_eval
TopoDS_Face.outer_wire = property(_outer_wire)
TopoDS_Face.make_halfspace = _make_halfspace
TopoDS_Face.area = property(_mass)

TopoDS_Shell.make_halfspace = _make_halfspace
TopoDS_Face.area = property(_mass)

TopoDS_Solid.outer_shell = property(_outer_shell)
TopoDS_Solid.volume = property(_mass)
