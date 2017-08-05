from OCC.BRep import BRep_Builder, BRep_Tool, BRep_Tool_Parameter
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.BRepClass3d import brepclass3d
from OCC.BRepTools import BRepTools_WireExplorer, breptools_OuterWire
from OCC.GeomAbs import (GeomAbs_BSplineCurve, GeomAbs_BSplineSurface,
                         GeomAbs_BezierCurve, GeomAbs_BezierSurface,
                         GeomAbs_Line, GeomAbs_Plane)
from OCC.GeomAdaptor import GeomAdaptor_Surface
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

from afem.geometry.entities import Line, Plane
from afem.geometry.methods.create import (create_nurbs_curve_from_occ,
                                          create_nurbs_surface_from_occ)

__all__ = ["ExploreShape", "ExploreWire", "ExploreFreeEdges"]


def _make_compound(shapes):
    """
    Make a compound from a list of shapes.
    """
    cp = TopoDS_Compound()
    builder = BRep_Builder()
    builder.MakeCompound(cp)
    for shape in shapes:
        builder.Add(cp, shape)
    return cp


class ExploreShape(object):
    """
    Explore shape.
    """

    @staticmethod
    def get_vertices(shape, as_compound=False):
        """
        Get vertices from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.

        :return: Vertices of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Vertex] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Vertex):
            if as_compound:
                return _make_compound([shape])
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_VERTEX)
        vertices = []
        while exp.More():
            vi = exp.Current()
            vertex = topods_Vertex(vi)
            vertices.append(vertex)
            exp.Next()
        if as_compound:
            return _make_compound(vertices)
        return vertices

    @staticmethod
    def get_edges(shape, as_compound=False, unique=True):
        """
        Get edges from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.
        :param bool unique: Option to return only unique edges.

        :return: Edges of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Edge] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Edge):
            if as_compound:
                return _make_compound([shape])
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
        if as_compound:
            return _make_compound(edges)
        return edges

    @staticmethod
    def get_wires(shape, as_compound=False):
        """
        Get wires from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.

        :return: Wires of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Wire] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Wire):
            if as_compound:
                return _make_compound([shape])
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_WIRE)
        wires = []
        while exp.More():
            wi = exp.Current()
            wire = topods_Wire(wi)
            wires.append(wire)
            exp.Next()
        if as_compound:
            return _make_compound(wires)
        return wires

    @staticmethod
    def get_faces(shape, as_compound=False):
        """
        Get faces from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.

        :return: Faces of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Face] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Face):
            if as_compound:
                return _make_compound([shape])
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_FACE)
        faces = []
        while exp.More():
            fi = exp.Current()
            face = topods_Face(fi)
            faces.append(face)
            exp.Next()
        if as_compound:
            return _make_compound(faces)
        return faces

    @staticmethod
    def get_shells(shape, as_compound=False):
        """
        Get shells from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.

        :return: Shells of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Shell] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Shell):
            if as_compound:
                return _make_compound([shape])
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_SHELL)
        shells = []
        while exp.More():
            si = exp.Current()
            shell = topods_Shell(si)
            shells.append(shell)
            exp.Next()
        if as_compound:
            return _make_compound(shells)
        return shells

    @staticmethod
    def get_solids(shape, as_compound=False):
        """
        Get solids from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.

        :return: Solids of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Solid] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Solid):
            if as_compound:
                return _make_compound([shape])
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_SOLID)
        solids = []
        while exp.More():
            si = exp.Current()
            solid = topods_Solid(si)
            solids.append(solid)
            exp.Next()
        if as_compound:
            return _make_compound(solids)
        return solids

    @staticmethod
    def get_compounds(shape, as_compound=False):
        """
        Get compounds from a shape.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param bool as_compound: Option to return as a compound rather than
            a list.

        :return: Compounds of shape.
        :rtype: list[OCC.TopoDS.TopoDS_Compound] or OCC.TopoDS.TopoDS_Compound
        """
        if isinstance(shape, TopoDS_Compound):
            if as_compound:
                return _make_compound([shape])
            return [shape]

        exp = TopExp_Explorer(shape, TopAbs_COMPOUND)
        compounds = []
        while exp.More():
            ci = exp.Current()
            compound = topods_Compound(ci)
            compounds.append(compound)
            exp.Next()
        if as_compound:
            return _make_compound(compounds)
        return compounds

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
    def curve_of_edge(edge):
        """
        Get the curve of the edge.

        :param OCC.TopoDS.TopoDS_Edge edge: The edge.

        :return: Underlying curve of edge.
        :rtype: afem.geometry.entities.Curve

        :raise RuntimeError: If an unsupported curve type is found.
        """
        # TODO Handle other curve types.
        adp_crv = BRepAdaptor_Curve(edge)
        if adp_crv.GetType() == GeomAbs_Line:
            gp_lin = adp_crv.Line()
            crv = Line(gp_lin)
            return crv
        elif adp_crv.GetType() in [GeomAbs_BezierCurve, GeomAbs_BSplineCurve]:
            crv = adp_crv.BSpline().GetObject()
            crv = create_nurbs_curve_from_occ(crv)
            return crv
        else:
            msg = 'Unsupported curve type.'
            raise RuntimeError(msg)

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
        occ_crv = geom_convert.BSplineCurve().GetObject()
        if not occ_crv:
            msg = 'Unsupported curve type.'
            raise RuntimeError(msg)
        crv = create_nurbs_curve_from_occ(occ_crv)
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

        :raise RuntimeError: If unsupported surface type is found.
        """
        hsrf = BRep_Tool.Surface(face)
        adp_srf = GeomAdaptor_Surface(hsrf)
        if adp_srf.GetType() == GeomAbs_Plane:
            gp_pln = adp_srf.Plane()
            return Plane(gp_pln)
        elif adp_srf.GetType() in [GeomAbs_BezierSurface,
                                   GeomAbs_BSplineSurface]:
            occ_srf = adp_srf.BSpline().GetObject()
            return create_nurbs_surface_from_occ(occ_srf)
        else:
            msg = 'Unsupported surface type.'
            raise RuntimeError(msg)

    @classmethod
    def surface_of_shape(cls, shape):
        """
        Get the surface of the largest face in the shape.

        :param shape:
        :return:
        """
        faces = cls.get_faces(shape)
        if not faces:
            return None

        # TODO Largest face method.
        f = cls.largest_face(faces)
        if not f:
            return None
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

    @classmethod
    def get_free_edges(cls, shape, as_compound=False):
        """
        Get the free edges of a shape.

        :param shape:
        :param bool as_compound:

        :return:
        """
        # TODO Use ExploreFreeEdges class.
        faces = ExploreShape.get_faces(shape)
        if not faces:
            return []

        compound = _make_compound(faces)
        tol = ExploreShape.get_tolerance(compound, 1)
        fb_tool = ShapeAnalysis_FreeBounds(compound, tol)
        closed_edges = ExploreShape.get_edges(fb_tool.GetClosedWires())
        open_edges = ExploreShape.get_edges(fb_tool.GetOpenWires())
        edges = closed_edges + open_edges
        if as_compound:
            return _make_compound(edges)
        return edges


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
        self._closed_wires = ExploreShape.get_wires(tool.GetClosedWires())
        self._open_wires = ExploreShape.get_wires(tool.GetOpenWires())
        self._edges = (ExploreShape.get_edges(self._closed_wires) +
                       ExploreShape.get_edges(self._open_wires))

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
