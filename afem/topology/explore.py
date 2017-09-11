from itertools import product

from OCC.BRep import BRep_Tool, BRep_Tool_Parameter
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.BRepBuilderAPI import BRepBuilderAPI_Copy
from OCC.BRepClass3d import brepclass3d
from OCC.BRepTools import BRepTools_WireExplorer, breptools_OuterWire
from OCC.GeomAbs import (GeomAbs_BSplineCurve, GeomAbs_BezierCurve)
from OCC.GeomConvert import GeomConvert_CompCurveToBSplineCurve
from OCC.ShapeAnalysis import (ShapeAnalysis_Edge, ShapeAnalysis_FreeBounds,
                               ShapeAnalysis_ShapeTolerance)
from OCC.TopAbs import (TopAbs_COMPOUND, TopAbs_EDGE, TopAbs_FACE,
                        TopAbs_SHELL, TopAbs_SOLID, TopAbs_VERTEX, TopAbs_WIRE)
from OCC.TopExp import TopExp_Explorer
from OCC.TopoDS import (TopoDS_Compound, TopoDS_Edge, TopoDS_Face,
                        TopoDS_Shell, TopoDS_Solid, TopoDS_Vertex,
                        TopoDS_Wire, topods_Compound, topods_Edge, topods_Face,
                        topods_Shell, topods_Solid, topods_Vertex, topods_Wire)

from afem.geometry.entities import Curve, NurbsCurve, Point, Surface
from afem.topology.props import AreaOfShapes

__all__ = ["ExploreShape", "ExploreWire", "ExploreFreeEdges"]


class ExploreShape(object):
    """
    Explore shape.
    """

    @staticmethod
    def get_vertices(shape, unique=True):
        """
        Get vertices from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool unique: Option to return only unique vertices.

        :return: Vertices of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Vertex]
        """
        # TODO Option for unique vertices.
        if isinstance(shape, TopoDS_Vertex):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_VERTEX)
        vertices = []
        while exp.More():
            vi = exp.Current()
            vertex = topods_Vertex(vi)
            if unique:
                is_unique = True
                for v in vertices:
                    if v.IsSame(vertex):
                        is_unique = False
                        break
                if is_unique:
                    vertices.append(vertex)
            else:
                vertices.append(vertex)
            exp.Next()
        return vertices

    @staticmethod
    def get_edges(shape, unique=True):
        """
        Get edges from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool unique: Option to return only unique edges.

        :return: Edges of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        if isinstance(shape, TopoDS_Edge):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_EDGE)
        edges = []
        while exp.More():
            ei = exp.Current()
            edge = topods_Edge(ei)
            if unique:
                is_unique = True
                for e in edges:
                    if e.IsSame(edge):
                        is_unique = False
                        break
                if is_unique:
                    edges.append(edge)
            else:
                edges.append(edge)
            exp.Next()
        return edges

    @staticmethod
    def get_wires(shape):
        """
        Get wires from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: Wires of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Wire]
        """
        if isinstance(shape, TopoDS_Wire):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_WIRE)
        wires = []
        while exp.More():
            wi = exp.Current()
            wire = topods_Wire(wi)
            wires.append(wire)
            exp.Next()
        return wires

    @staticmethod
    def get_faces(shape):
        """
        Get faces from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: Faces of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Face]
        """
        if isinstance(shape, TopoDS_Face):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_FACE)
        faces = []
        while exp.More():
            fi = exp.Current()
            face = topods_Face(fi)
            faces.append(face)
            exp.Next()
        return faces

    @staticmethod
    def get_shells(shape):
        """
        Get shells from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: Shells of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Shell]
        """
        if isinstance(shape, TopoDS_Shell):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_SHELL)
        shells = []
        while exp.More():
            si = exp.Current()
            shell = topods_Shell(si)
            shells.append(shell)
            exp.Next()
        return shells

    @staticmethod
    def get_solids(shape):
        """
        Get solids from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: Solids of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Solid]
        """
        if isinstance(shape, TopoDS_Solid):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_SOLID)
        solids = []
        while exp.More():
            si = exp.Current()
            solid = topods_Solid(si)
            solids.append(solid)
            exp.Next()
        return solids

    @staticmethod
    def get_compounds(shape):
        """
        Get compounds from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: Compounds of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Compound]
        """
        if isinstance(shape, TopoDS_Compound):
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_COMPOUND)
        compounds = []
        while exp.More():
            ci = exp.Current()
            compound = topods_Compound(ci)
            compounds.append(compound)
            exp.Next()
        return compounds

    @staticmethod
    def get_shared_edges(shape1, shape2):
        """
        Get the shared edges between the two shapes.

        :param TopoDS.TopoDS_Shape shape1: The first shape.
        :param TopoDS.TopoDS_Shape shape2: The second shape.

        :rtype: List of shared edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        edges1 = ExploreShape.get_edges(shape1)
        edges2 = ExploreShape.get_edges(shape2)
        if not edges1 or not edges2:
            return []

        shared_edges = []
        for e1, e2 in product(edges1, edges2):
            if e1.IsSame(e2):
                unique = True
                for ei in shared_edges:
                    if ei.IsSame(e1):
                        unique = False
                        break
                if unique:
                    shared_edges.append(e1)

        return shared_edges

    @staticmethod
    def outer_wire(face):
        """
        Get outer wire of face.

        :param OCC.TopoDS.TopoDS_Face face: The face.

        :return: Outer wire.
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return breptools_OuterWire(face)

    @staticmethod
    def outer_shell(solid):
        """
        Get the outer shell of the solid.

        :param OCC.TopoDS.TopoDS_Solid solid: The solid.

        :return: Outer shell.
        :rtype: OCC.TopoDS.TopoDS_Shell
        """
        return brepclass3d.OuterShell(solid)

    @staticmethod
    def pnt_of_vertex(vertex):
        """
        Get the underlying point of the vertex.

        :param OCC.TopoDS.TopoDS_Vertex vertex: The vertex.

        :return: The point.
        :rtype: afem.geometry.entities.Point
        """
        gp_pnt = BRep_Tool.Pnt(vertex)
        return Point(gp_pnt.X(), gp_pnt.Y(), gp_pnt.Z())

    @staticmethod
    def curve_of_edge(edge):
        """
        Get the curve of the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: Underlying curve of edge.
        :rtype: afem.geometry.entities.Curve
        """
        h_crv = BRep_Tool.Curve(edge)
        return Curve(h_crv[0])

    @staticmethod
    def curve_of_wire(wire):
        """
        Get the curve of the wire. The edges are concatenated so the
        resulting curve may be C0 continuous.

        :param OCC.TopoDS.TopoDS_Wire wire: The wire.

        :return: Concatenated curve of wire.
        :rtype: afem.geometry.entities.Curve

        :raise RuntimeError: If an unsupported curve type is found.
        """
        geom_convert = GeomConvert_CompCurveToBSplineCurve()
        exp = BRepTools_WireExplorer(wire)
        while exp.More():
            e = topods_Edge(exp.Current())
            exp.Next()
            adp_crv = BRepAdaptor_Curve(e)
            tol = ExploreShape.get_tolerance(e, 1)
            # TODO Handle curves other than BSpline
            if adp_crv.GetType() not in [GeomAbs_BSplineCurve,
                                         GeomAbs_BezierCurve]:
                continue
            geom_convert.Add(adp_crv.BSpline(), tol)
        h_crv = geom_convert.BSplineCurve()
        if h_crv.IsNull():
            msg = 'Unsupported curve type.'
            raise RuntimeError(msg)
        crv = NurbsCurve(h_crv)
        return crv

    @classmethod
    def curve_of_shape(cls, shape):
        """
        Get the curve of the shape if possible.

        :param shape: The edge or wire.
        :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire

        :return: The underlying curve of the edge or wire.
        :rtype: afem.geometry.entities.Curve

        :raise TypeError: If *shape* is not an edge or wire.
        """
        if shape.ShapeType() == TopAbs_EDGE:
            return cls.curve_of_edge(shape)
        elif shape.ShapeType() == TopAbs_WIRE:
            return cls.curve_of_wire(shape)
        else:
            msg = 'Invalid shape.'
            raise TypeError(msg)

    @staticmethod
    def surface_of_face(face):
        """
        Get the surface of the face.

        :param OCC.TopoDS.TopoDS_Face face: The face.

        :return: Underlying surface of face.
        :rtype: afem.geometry.entities.Surface
        """
        hsrf = BRep_Tool.Surface(face)
        return Surface(hsrf)

    @classmethod
    def surface_of_shape(cls, shape):
        """
        Get the surface of the largest face in the shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: The surface.
        :rtype: afem.geometry.entities.Surface
        """
        faces = cls.get_faces(shape)
        if not faces:
            return None

        f = AreaOfShapes(faces).largest_shape
        return cls.surface_of_face(f)

    @staticmethod
    def first_vertex(edge):
        """
        Return the first vertex of the edge considering orientation.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: The first vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return ShapeAnalysis_Edge().FirstVertex(edge)

    @staticmethod
    def last_vertex(edge):
        """
        Return the last vertex of the edge considering orientation.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: The last vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return ShapeAnalysis_Edge().LastVertex(edge)

    @classmethod
    def vertices(cls, edge):
        """
        Return the first and last vertex of the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: The first and last vertices (v1, v2).
        :rtype: tuple(OCC.TopoDS.TopoDS_Vertex)
        """
        return cls.first_vertex(edge), cls.last_vertex(edge)

    @staticmethod
    def parameter(vertex, edge, face=None):
        """
        Return the parameter of the vertex on the edge.

        :param OCC.TopoDS.TopoDS_Vertex vertex: The vertex.
        :param OCC.TopoDS.TopoDS_Edge edge: The edge.
        :param OCC.TopoDS.TopoDS_Face face: The face.

        :return: The parameter.
        :rtype: float
        """
        if not face:
            return BRep_Tool_Parameter(vertex, edge)
        else:
            return BRep_Tool_Parameter(vertex, edge, face)

    @classmethod
    def parameters(cls, edge, face=None):
        """
        Return the first and last parameters on the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.
        :param OCC.TopoDS.TopoDS_Face face: The face.

        :return: The parameters (u1, u2).
        :rtype: tuple(float)
        """
        v1, v2 = cls.vertices(edge)
        u1 = cls.parameter(v1, edge, face)
        u2 = cls.parameter(v2, edge, face)
        return u1, u2

    @staticmethod
    def get_tolerance(shape, mode=0):
        """
        Compute the global tolerance of the shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param int mode: Average (0), maximal (1), minimal (2)

        :return: The tolerance.
        :rtype: float
        """
        tol = ShapeAnalysis_ShapeTolerance()
        tol.AddTolerance(shape)
        return tol.GlobalTolerance(mode)

    @staticmethod
    def copy_shape(shape, copy_geom=True):
        """
        Copy a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape,
        :param bool copy_geom: Option to copy geometry.

        :return: The copied shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return BRepBuilderAPI_Copy(shape, copy_geom).Shape()


class ExploreWire(object):
    """
    Explore the edges of a wire.

    :param OCC.TopoDS.TopoDS_Wire wire: The wire.
    :param OCC.TopoDS.TopoDS_Face face: The face.

    Usage:

    >>> from afem.topology import ExploreWire, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> explorer = ExploreWire(wire)
    >>> explorer.nedges
    4
    """

    def __init__(self, wire, face=None):
        if face is None:
            explorer = BRepTools_WireExplorer(wire)
        else:
            explorer = BRepTools_WireExplorer(wire, face)

        edges = []
        while explorer.More():
            ei = topods_Edge(explorer.Current())
            edges.append(ei)
            explorer.Next()
        self._edges = edges

    @property
    def nedges(self):
        """
        :return: Number of edges.
        :rtype: int
        """
        return len(self._edges)

    @property
    def edges(self):
        """
        :return: The ordered edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return self._edges


class ExploreFreeEdges(object):
    """
    Explore the free bounds of a shape.

    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    """

    def __init__(self, shape):
        tool = ShapeAnalysis_FreeBounds(shape)
        cmp_closed_wires = tool.GetClosedWires()
        cmp_open_wires = tool.GetOpenWires()
        self._closed_wires = ExploreShape.get_wires(cmp_closed_wires)
        self._open_wires = ExploreShape.get_wires(cmp_open_wires)
        self._edges = (ExploreShape.get_edges(cmp_closed_wires) +
                       ExploreShape.get_edges(cmp_open_wires))

    @property
    def closed_wires(self):
        """
        :return: Closed wires of free edges.
        :rtype: list[OCC.TopoDS.TopoDS_Wire]
        """
        return self._closed_wires

    @property
    def open_wires(self):
        """
        :return: Open wires of free edges.
        :rtype: list[OCC.TopoDS.TopoDS_Wire]
        """
        return self._open_wires

    @property
    def free_edges(self):
        """
        :return: All free edges of the shape.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return self._edges


if __name__ == "__main__":
    import doctest

    doctest.testmod()
