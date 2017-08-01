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
from .check import CheckShape
from ..geometry.methods.create import (create_nurbs_curve_from_occ,
                                       create_nurbs_surface_from_occ)


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

        :return:
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

        :param shape:
        :param bool unique: Unique edges based on TShape and Location.
        :param as_compound:

        :return:
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

        :return:
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

        :return:
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

        :return:
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

        :return:
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

        :return:
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

        :param face:

        :return:
        """
        face = CheckShape.to_face(face)
        if not face:
            return None

        return breptools_OuterWire(face)

    @staticmethod
    def outer_shell(solid):
        """
        Get the outer shell of the solid.

        :param solid:

        :return:
        """
        solid = CheckShape.to_solid(solid)
        if not solid:
            return None

        return brepclass3d.OuterShell(solid)

    @staticmethod
    def curve_of_edge(edge):
        """
        Get the curve of the edge.

        :param edge:

        :return:
        """
        edge = CheckShape.to_edge(edge)
        if not edge:
            return None

        adp_crv = BRepAdaptor_Curve(edge)
        if adp_crv.GetType() == GeomAbs_Line:
            gp_lin = adp_crv.Line()
            crv = Line(gp_lin)
            return crv
        # TODO Handle other curve types.
        crv = adp_crv.BSpline().GetObject()
        crv = create_nurbs_curve_from_occ(crv)
        return crv

    @staticmethod
    def curve_of_wire(wire):
        """
        Get the curve of the wire.

        :return:
        """
        wire = CheckShape.to_wire(wire)
        if not wire:
            return None

        geom_convert = GeomConvert_CompCurveToBSplineCurve()
        exp = ExploreShape.wire_explorer(wire)
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
            return None
        crv = create_nurbs_curve_from_occ(occ_crv)
        return crv

    @classmethod
    def curve_of_shape(cls, shape):
        """
        Get the curve of the shape if possible.

        :param shape:
        :return:
        """
        shape = CheckShape.to_shape(shape)
        if isinstance(shape, TopoDS_Edge):
            return cls.curve_of_edge(shape)
        elif isinstance(shape, TopoDS_Wire):
            return cls.curve_of_wire(shape)
        else:
            return None

    @staticmethod
    def surface_of_face(face):
        """
        Get the surface of the face.

        :param face:

        :return:
        """
        face = CheckShape.to_face(face)
        if not face:
            return None

        hsrf = BRep_Tool.Surface(face)
        adp_srf = GeomAdaptor_Surface(hsrf)
        if adp_srf.GetType() == GeomAbs_Plane:
            gp_pln = adp_srf.Plane()
            return Plane(gp_pln)
        elif adp_srf.GetType() in [GeomAbs_BezierSurface,
                                   GeomAbs_BSplineSurface]:
            occ_srf = adp_srf.BSpline().GetObject()
            return create_nurbs_surface_from_occ(occ_srf)

        return None

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
    def wire_explorer(wire, face=None):
        """
        Create a wire explorer.

        :param wire:
        :param face:

        :return:
        """
        wire = CheckShape.to_wire(wire)
        face = CheckShape.to_face(face)
        if not wire:
            return None

        # TODO Use custom ExploreWire class
        if not face:
            return BRepTools_WireExplorer(wire)
        else:
            return BRepTools_WireExplorer(wire, face)

    @staticmethod
    def first_vertex(edge):
        """
        Return the first vertex of the edge considering orientation.

        :param edge:

        :return:
        """
        edge = CheckShape.to_edge(edge)
        if not edge:
            return None

        return ShapeAnalysis_Edge().FirstVertex(edge)

    @staticmethod
    def last_vertex(edge):
        """
        Return the last vertex of the edge considering orientation.

        :param edge:

        :return:
        """
        edge = CheckShape.to_edge(edge)
        if not edge:
            return None

        return ShapeAnalysis_Edge().LastVertex(edge)

    @classmethod
    def vertices(cls, edge):
        """
        Return the first and last vertex of the edge.

        :param edge:

        :return:
        """
        edge = CheckShape.to_edge(edge)
        if not edge:
            return None, None

        return cls.first_vertex(edge), cls.last_vertex(edge)

    @staticmethod
    def parameter(vertex, edge, face=None):
        """
        Return the parameter of the vertex on the edge.

        :param vertex:
        :param edge:
        :param face:

        :return:
        """
        vertex = CheckShape.to_vertex(vertex)
        edge = CheckShape.to_edge(edge)
        face = CheckShape.to_face(face)
        if not vertex or not edge:
            return None
        if not face:
            return BRep_Tool_Parameter(vertex, edge)
        else:
            return BRep_Tool_Parameter(vertex, edge, face)

    @classmethod
    def parameters(cls, edge, face=None):
        """
        Return the first and last parameters on the edge.

        :param edge:
        :param face:

        :return:
        """
        v1, v2 = cls.vertices(edge)
        u1 = cls.parameter(v1, edge, face)
        u2 = cls.parameter(v2, edge, face)
        return u1, u2

    @staticmethod
    def get_tolerance(shape, mode=0):
        """
        Compute the global tolerance of the shape.

        :param shape:
        :param mode: Average (0), maximal (1), minimal (2)

        :return:
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
    pass


class ExploreFreeEdges(object):
    pass
