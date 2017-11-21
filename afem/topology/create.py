from OCCT.BRep import BRep_Builder, BRep_Tool
from OCCT.BRepAdaptor import BRepAdaptor_CompCurve, BRepAdaptor_Curve
from OCCT.BRepAlgo import BRepAlgo
from OCCT.BRepAlgoAPI import BRepAlgoAPI_Section
from OCCT.BRepBuilderAPI import (BRepBuilderAPI_FindPlane,
                                 BRepBuilderAPI_MakeEdge,
                                 BRepBuilderAPI_MakeFace,
                                 BRepBuilderAPI_MakePolygon,
                                 BRepBuilderAPI_MakeShell,
                                 BRepBuilderAPI_MakeSolid,
                                 BRepBuilderAPI_MakeVertex,
                                 BRepBuilderAPI_MakeWire, BRepBuilderAPI_Sewing)
from OCCT.BRepMesh import BRepMesh_IncrementalMesh
from OCCT.BRepOffsetAPI import (BRepOffsetAPI_MakeOffset)
from OCCT.BRepPrimAPI import (BRepPrimAPI_MakeCylinder,
                              BRepPrimAPI_MakeHalfSpace, BRepPrimAPI_MakePrism)
from OCCT.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
# from OCCT.GEOMAlgo import GEOMAlgo_Splitter
from OCCT.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCCT.GeomAbs import GeomAbs_Arc, GeomAbs_Intersection, GeomAbs_Tangent
from OCCT.ShapeAnalysis import ShapeAnalysis_FreeBounds
from OCCT.ShapeBuild import ShapeBuild_ReShape
from OCCT.TopAbs import (TopAbs_COMPOUND, TopAbs_COMPSOLID, TopAbs_EDGE,
                         TopAbs_FACE, TopAbs_SHELL, TopAbs_SOLID, TopAbs_WIRE)
from OCCT.TopTools import TopTools_HSequenceOfShape
from OCCT.TopoDS import (TopoDS_Compound, TopoDS_Face, TopoDS_Shape,
                         TopoDS_Shell, TopoDS)
from numpy import ceil

from afem.geometry.check import CheckGeom
from afem.geometry.create import PlaneByApprox
from afem.geometry.entities import Plane, Point
from afem.topology.bop import IntersectShapes
from afem.topology.check import CheckShape
from afem.topology.explore import ExploreShape

__all__ = ["VertexByPoint", "EdgeByPoints", "EdgeByVertices", "EdgeByCurve",
           "EdgeByDrag", "EdgeByWireConcat", "WireByEdges",
           "WiresByConnectedEdges",
           "WireByPlanarOffset", "WiresByShape", "WireByPoints", "WireBySplit",
           "WireByConcat", "FaceBySurface",
           "FaceByPlane", "FaceByPlanarWire", "FaceByDrag", "ShellBySurface",
           "ShellByFaces", "ShellBySewing",
           "SolidByShell", "ShellByDrag", "SolidByPlane", "SolidByDrag",
           "SolidByCylinder",
           "CompoundByShapes", "HalfspaceByShape", "HalfspaceBySurface",
           "ShapeByFaces", "ShapeByDrag", "PointAlongShape",
           "PointsAlongShapeByNumber", "PointsAlongShapeByDistance",
           "PlaneByEdges", "PlaneByIntersectingShapes"]

_occ_join = {'a': GeomAbs_Arc,
             'arc': GeomAbs_Arc,
             't': GeomAbs_Tangent,
             'tangent': GeomAbs_Tangent,
             'i': GeomAbs_Intersection,
             'intersection': GeomAbs_Intersection}


# VERTEX ----------------------------------------------------------------------

class VertexByPoint(object):
    """
    Create a vertex using a point.

    :param point_like pnt: The point.

    Usage:

    >>> from afem.topology import VertexByPoint
    >>> builder = VertexByPoint((1., 2., 3.))
    >>> v = builder.vertex
    """

    def __init__(self, pnt):
        # Build
        pnt = CheckGeom.to_point(pnt)
        builder = BRepBuilderAPI_MakeVertex(pnt)
        v = builder.Vertex()
        self._v = v

    @property
    def vertex(self):
        """
        :return: The vertex.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v


# EDGE ------------------------------------------------------------------------

class EdgeByPoints(object):
    """
    Create an edge between two points.

    :param point_like p1: The first point.
    :param point_like p2: The second point.

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> builder = EdgeByPoints((0., 0., 0.), (10., 0., 0.))
    >>> e = builder.edge
    >>> v1 = builder.vertex1
    >>> v2 = builder.vertex2
    """

    def __init__(self, p1, p2):
        # Build
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        builder = BRepBuilderAPI_MakeEdge(p1, p2)
        self._e = builder.Edge()
        self._v1 = builder.Vertex1()
        self._v2 = builder.Vertex2()

    @property
    def edge(self):
        """
        :return: The edge.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._e

    @property
    def vertex1(self):
        """
        :return: The first vertex.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v1

    @property
    def vertex2(self):
        """
        :return: The second vertex.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v2


class EdgeByVertices(object):
    """
    Create an edge between two vertices.

    :param OCCT.TopoDS.TopoDS_Vertex v1: The first vertex.
    :param OCCT.TopoDS.TopoDS_Vertex v2: The second vertex.

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> builder = EdgeByPoints((0., 0., 0.), (10., 0., 0.))
    >>> e = builder.edge
    """

    def __init__(self, v1, v2):
        # Build
        e = BRepBuilderAPI_MakeEdge(v1, v2).Edge()
        self._e = e

    @property
    def edge(self):
        """
        :return: The edge.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._e


class EdgeByCurve(object):
    """
    Create an edge using a curve.

    :param afem.geometry.entities.Curve crv: The curve.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.topology import EdgeByCurve
    >>> c = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> builder = EdgeByCurve(c)
    >>> e = builder.edge
    >>> v1 = builder.vertex1
    >>> v2 = builder.vertex2
    """

    def __init__(self, crv):
        # Build
        builder = BRepBuilderAPI_MakeEdge(crv.handle)
        self._e = builder.Edge()
        self._v1 = builder.Vertex1()
        self._v2 = builder.Vertex2()

    @property
    def edge(self):
        """
        :return: The edge.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._e

    @property
    def vertex1(self):
        """
        :return: The first vertex.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v1

    @property
    def vertex2(self):
        """
        :return: The second vertex.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v2


class EdgeByDrag(object):
    """
    Create an edge by dragging a vertex along a vector.

    :param OCCT.TopoDS.TopoDS_Vertex vertex: The vertex.
    :param vector_like v: The vector to drag the shape.

    Usage:

    >>> from afem.topology import EdgeByDrag, VertexByPoint
    >>> vertex = VertexByPoint((0., 0., 0.)).vertex
    >>> builder = EdgeByDrag(vertex, (1., 0., 0.))
    >>> e = builder.edge
    >>> v1 = builder.first_vertex
    >>> v2 = builder.last_vertex
    """

    def __init__(self, vertex, v):
        v = CheckGeom.to_vector(v)
        builder = BRepPrimAPI_MakePrism(vertex, v)
        self._e = TopoDS.Edge_(builder.Shape())
        self._v1 = TopoDS.Vertex_(builder.FirstShape())
        self._v2 = TopoDS.Vertex_(builder.LastShape())

    @property
    def edge(self):
        """
        :return: The edge.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._e

    @property
    def first_vertex(self):
        """
        :return: The vertex at the bottom of the edge.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v1

    @property
    def last_vertex(self):
        """
        :return: The vertex at the top of the edge.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._v2


class EdgeByWireConcat(object):
    """
    Create an edge by concatenating all the edges of a wire. The edge may
    have C0 continuity.

    :param OCCT.TopoDS.TopoDS_Wire wire: The wire.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.topology import EdgeByCurve, WireByEdges, EdgeByWireConcat
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(10., 0., 0.), (11., 1., 0.)]).curve
    >>> e1 =  EdgeByCurve(c1).edge
    >>> e2 =  EdgeByCurve(c2).edge
    >>> w = WireByEdges(e1, e2).wire
    >>> edge = EdgeByWireConcat(w).edge
    """

    def __init__(self, wire):
        self._edge = BRepAlgo.ConcatenateWireC0_(wire)

    @property
    def edge(self):
        """
        :return: The edge.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._edge


# WIRE ------------------------------------------------------------------------


class WireByEdges(object):
    """
    Create a wire using up to four edges.

    :param OCCT.TopoDS.TopoDS_Edge e1: The first edge.
    :param OCCT.TopoDS.TopoDS_Edge e2: The second edge.
    :param OCCT.TopoDS.TopoDS_Edge e3: The third edge.
    :param OCCT.TopoDS.TopoDS_Edge e4: The fourth edge.

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.topology import EdgeByCurve, WireByEdges
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(10., 0., 0.), (11., 0., 0.)]).curve
    >>> e1 =  EdgeByCurve(c1).edge
    >>> e2 =  EdgeByCurve(c2).edge
    >>> builder = WireByEdges(e1, e2)
    >>> w = builder.wire
    >>> e = builder.last_edge
    >>> v = builder.last_vertex
    """

    def __init__(self, e1, e2=None, e3=None, e4=None):
        # Build
        builder = BRepBuilderAPI_MakeWire(e1)
        for e in [e2, e3, e4]:
            if e is not None and e.ShapeType() == TopAbs_EDGE:
                builder.Add(e)

        self._w = builder.Wire()
        self._last_e = builder.Edge()
        self._last_v = builder.Vertex()

    @property
    def wire(self):
        """
        :return: The wire.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._w

    @property
    def last_edge(self):
        """
        :return: The last edge added to the wire.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._last_e

    @property
    def last_vertex(self):
        """
        :return: The last vertex added to the wire.
        :rtype: OCCT.TopoDS.TopoDS_Vertex
        """
        return self._last_v


class WiresByConnectedEdges(object):
    """
    Create wires from a list of unsorted edges.

    :param list[OCCT.TopoDS.TopoDS_Edge] edges: List of edges.
    :param float tol: Connection tolerance. If *None* if provided then the
        maximum tolerance of all edge will be used.
    :param bool shared: Option to use only shared vertices to connect edges.
        If *False* then geometric coincidence will be also checked.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.topology import EdgeByCurve, WiresByConnectedEdges
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(10., 0., 0.), (11., 0., 0.)]).curve
    >>> e1 =  EdgeByCurve(c1).edge
    >>> e2 =  EdgeByCurve(c2).edge
    >>> builder = WiresByConnectedEdges([e1, e2])
    >>> builder.nwires
    1
    """

    def __init__(self, edges, tol=None, shared=False):
        # Build
        hedges = TopTools_HSequenceOfShape()
        for e in edges:
            hedges.Append(e)

        if tol is None:
            tol = max([ExploreShape.global_tolerance(e, 1) for e in edges])

        hwires = ShapeAnalysis_FreeBounds.ConnectEdgesToWires_(hedges, tol,
                                                               shared)

        wires = []
        for i in range(1, hwires.Length() + 1):
            w = TopoDS.Wire_(hwires.Value(i))
            wires.append(w)

        self._wires = wires

    @property
    def nwires(self):
        """
        :return: Number of wires.
        :rtype: int
        """
        return len(self._wires)

    @property
    def wires(self):
        """
        :return: The wires.
        :rtype: list[OCCT.TopoDS.TopoDS_Wire]
        """
        return self._wires


class WireByPlanarOffset(object):
    """
    Create a wire by offsetting a planar wire or face.

    :param spine: The wire to offset. If a face is provided the outer wire
        will be used.
    :type spine: OCCT.TopoDS.TopoDS_Wire or OCCT.TopoDS.TopoDS_Face
    :param float distance: Offset distance in the plane.
    :param float altitude: Offset altitude normal to the plane.
    :param str join: Join type ('arc', 'tangent', 'intersection').
    """

    def __init__(self, spine, distance, altitude=0., join='arc',
                 is_open=False):
        try:
            join = join.lower()
            join = _occ_join[join]
        except (KeyError, AttributeError):
            join = GeomAbs_Arc

        offset = BRepOffsetAPI_MakeOffset(spine, join, is_open)
        offset.Perform(distance, altitude)
        shp = offset.Shape()
        w = CheckShape.to_wire(shp)
        self._w = w

    @property
    def wire(self):
        """
        :return: The wire.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._w


class WiresByShape(WiresByConnectedEdges):
    """
    Create wires by connecting all the edges of a shape. This method gathers
    all the unique edges of a shape and then uses
    :class:`.WiresByConnectedEdges`.

    :param OCCT.TopoDS.TopoDS_Shape: The shape.

    :raise ValueError: If no edges are found in the shape.
    """

    def __init__(self, shape):
        edges = ExploreShape.get_edges(shape)
        if not edges:
            msg = 'There are no edges in the shape to connect.'
            raise ValueError(msg)

        super(WiresByShape, self).__init__(edges)


class WireByPoints(object):
    """
    Create a polygonal wire by connecting points.

    :param list[point_like] pnts: The ordered points.
    :param bool close: Option to close the wire.

    Usage:

    >>> from afem.topology import WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> builder = WireByPoints([p1, p2, p3, p4], True)
    >>> wire = builder.wire
    >>> wire.Closed()
    True
    """

    def __init__(self, pnts, close=False):
        builder = BRepBuilderAPI_MakePolygon()
        for p in pnts:
            p = CheckGeom.to_point(p)
            builder.Add(p)
        if close:
            builder.Close()
        self._w = builder.Wire()

    @property
    def wire(self):
        """
        :return: The wire.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._w


class WireBySplit(object):
    """
    Split a wire with a shape.

    :param OCCT.TopoDS.TopoDS_Wire wire: The wire.
    :param OCCT.TopoDS.TopoDS_Shape splitter: The splitter shape.
    """

    def __init__(self, wire, splitter):
        # Split algorithm
        bop = GEOMAlgo_Splitter()
        bop.AddArgument(wire)
        bop.AddTool(splitter)
        bop.Perform()
        if bop.ErrorStatus() != 0:
            msg = 'Failed to split wire.'
            raise RuntimeError(msg)

        # Replace edges in wire
        reshape = ShapeBuild_ReShape()
        performed = False
        for old_edge in ExploreShape.get_edges(wire):
            # Check deleted
            if bop.IsDeleted(old_edge):
                reshape.Remove(old_edge)
                performed = True
                break

            # Check modified
            modified = bop.Modified(old_edge)
            if modified.IsEmpty():
                continue

            # Put modified edges into compound
            new_edges = []
            while not modified.IsEmpty():
                shape = modified.First()
                new_edges.append(shape)
                modified.RemoveFirst()
            if not new_edges:
                continue

            # Replace old edge with new.
            new_edges = CompoundByShapes(new_edges).compound
            reshape.Replace(old_edge, new_edges)
            performed = True

        # Apply substitution.
        if not performed:
            self._wire = wire
        else:
            new_wire = reshape.Apply(wire)
            self._wire = TopoDS.Wire_(new_wire)

    @property
    def wire(self):
        """
        :return: The split wire.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._wire


class WireByConcat(object):
    """
    Create a wire by concatenating all the edges of the wire. The wire may
    have C0 continuity.

    :param OCCT.TopoDS.TopoDS_Wire wire: The wire.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.topology import EdgeByCurve, WireByEdges, WireByConcat
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(10., 0., 0.), (11., 1., 0.)]).curve
    >>> e1 =  EdgeByCurve(c1).edge
    >>> e2 =  EdgeByCurve(c2).edge
    >>> w = WireByEdges(e1, e2).wire
    >>> wire = WireByConcat(w).wire
    """

    def __init__(self, wire):
        edge = BRepAlgo.ConcatenateWireC0_(wire)
        self._wire = BRepBuilderAPI_MakeWire(edge).Wire()

    @property
    def wire(self):
        """
        :return: The wire.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._wire


# FACE ------------------------------------------------------------------------

class FaceBySurface(object):
    """
    Create a face from a surface.

    :param afem.geometry.entities.Surface srf: The surface.
    :param float tol: Tolerance for resolution of degenerate edges.

    Usage:

    >>> from afem.geometry import PlaneByNormal
    >>> from afem.topology import FaceBySurface
    >>> pln = PlaneByNormal().plane
    >>> f = FaceBySurface(pln).face
    """

    def __init__(self, srf, tol=1.0e-7):
        builder = BRepBuilderAPI_MakeFace(srf.handle, tol)
        self._f = builder.Face()

    @property
    def face(self):
        """
        :return: The face.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._f


class FaceByPlane(object):
    """
    Create a finite face from a plane.

    :param afem.geometry.entities.Plane pln: The plane.
    :param float umin: Minimum u-parameter.
    :param float umax: Maximum u-parameter.
    :param float vmin: Minimum v-parameter.
    :param float vmax: Maximum v-parameter.

    Usage:

    >>> from afem.geometry import PlaneByNormal
    >>> from afem.topology import FaceByPlane
    >>> pln = PlaneByNormal().plane
    >>> f = FaceByPlane(pln, -1., 1., -1., 1.).face
    """

    def __init__(self, pln, umin, umax, vmin, vmax):
        builder = BRepBuilderAPI_MakeFace(pln.handle.Pln(), umin, umax, vmin,
                                          vmax)
        self._f = builder.Face()

    @property
    def face(self):
        """
        :return: The face.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._f


class FaceByPlanarWire(object):
    """
    Create a face from a planar wire.

    :param OCCT.TopoDS.TopoDS_Wire wire: The wire.

    Usage:

    >>> from afem.topology import FaceByPlanarWire, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> builder = FaceByPlanarWire(wire)
    >>> face = builder.face
    """

    def __init__(self, wire):
        self._f = BRepBuilderAPI_MakeFace(wire, True).Face()

    @property
    def face(self):
        """
        :return: The face.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._f


class FaceByDrag(object):
    """
    Create a face by dragging an edge along a vector.

    :param OCCT.TopoDS.TopoDS_Edge edge: The edge.
    :param vector_like v: The vector to drag the shape.

    Usage:

    >>> from afem.topology import EdgeByDrag, FaceByDrag, VertexByPoint
    >>> vertex = VertexByPoint((0., 0., 0.)).vertex
    >>> e = EdgeByDrag(vertex, (1., 0., 0.)).edge
    >>> builder = FaceByDrag(e, (0., 1., 0.))
    >>> f = builder.face
    >>> e1 = builder.first_edge
    >>> e2 = builder.last_edge
    """

    def __init__(self, edge, v):
        v = CheckGeom.to_vector(v)
        builder = BRepPrimAPI_MakePrism(edge, v)
        self._f = TopoDS.Face_(builder.Shape())
        self._e1 = TopoDS.Edge_(builder.FirstShape())
        self._e2 = TopoDS.Edge_(builder.LastShape())

    @property
    def face(self):
        """
        :return: The face.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._f

    @property
    def first_edge(self):
        """
        :return: The edge at the bottom of the face.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._e1

    @property
    def last_edge(self):
        """
        :return: The edge at the top of the face.
        :rtype: OCCT.TopoDS.TopoDS_Edge
        """
        return self._e2


# SHELL -----------------------------------------------------------------------

class ShellBySurface(object):
    """
    Create a shell from a surface.

    :param afem.geometry.entities.Surface srf: The surface.

    Usage:

    >>> from afem.geometry import PlaneByNormal
    >>> from afem.topology import ShellBySurface
    >>> pln = PlaneByNormal().plane
    >>> shell = ShellBySurface(pln).shell
    """

    def __init__(self, srf):
        builder = BRepBuilderAPI_MakeShell(srf.handle)
        self._shell = builder.Shell()

    @property
    def shell(self):
        """
        :return: The shell.
        :rtype: OCCT.TopoDS.TopoDS_Shell
        """
        return self._shell


class ShellByFaces(object):
    """
    Create a shell from connected faces. This method initializes a shell an
    then simply adds the faces to it. The faces should already have shared
    edges. This is not checked.

    :param list[OCCT.TopoDS.TopoDS_Face] faces: List of faces.

    Usage:

    >>> from afem.topology import FaceByPlanarWire, ShellByFaces, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> face = FaceByPlanarWire(wire).face
    >>> shell = ShellByFaces([face]).shell
    """

    def __init__(self, faces):
        shell = TopoDS_Shell()
        builder = BRep_Builder()
        builder.MakeShell(shell)
        for face in faces:
            builder.Add(shell, face)
        self._shell = shell

    @property
    def shell(self):
        """
        :return: The shell.
        :rtype: OCCT.TopoDS.TopoDS_Shell
        """
        return self._shell


class ShellByDrag(object):
    """
    Create a shell by dragging a wire along a vector.

    :param OCCT.TopoDS.TopoDS_Wire wire: The wire.
    :param vector_like v: The vector to drag the shape.

    Usage:

    >>> from afem.topology import ShellByDrag, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> builder = ShellByDrag(wire, (0., 0., 1.))
    >>> shell = builder.shell
    >>> w1 = builder.first_wire
    >>> w2 = builder.last_wire
    """

    def __init__(self, wire, v):
        v = CheckGeom.to_vector(v)
        builder = BRepPrimAPI_MakePrism(wire, v)
        self._shell = TopoDS.Shell_(builder.Shape())
        self._w1 = TopoDS.Wire_(builder.FirstShape())
        self._w2 = TopoDS.Wire_(builder.LastShape())

    @property
    def shell(self):
        """
        :return: The shell.
        :rtype: OCCT.TopoDS.TopoDS_Shell
        """
        return self._shell

    @property
    def first_wire(self):
        """
        :return: The wire at the bottom of the shell.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._w1

    @property
    def last_wire(self):
        """
        :return: The wire at the top of the shell.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._w2


class ShellBySewing(object):
    """
    Create a shell by sewing faces.

    :param list[OCCT.TopoDS.TopoDS_Face] faces: The faces.
    :param float tol: Sewing tolerance. If *None* the maximum tolerance of all
        the faces is used.
    :param bool cut_free_edges: Option for cutting of free edges.
    :param bool non_manifold: Option for non-manifold processing.

    :raise RuntimeError: If an invalid shape type results.

    Usage:

    >>> from afem.topology import *
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> w1 = WireByPoints([p1, p2, p3, p4], True).wire
    >>> f1 = FaceByPlanarWire(w1).face
    >>> p1 = (0., 0., 0.)
    >>> p2 = (0., 0., 1.)
    >>> p3 = (0., 1., 1.)
    >>> p4 = (0., 1., 0.)
    >>> w2 = WireByPoints([p1, p2, p3, p4], True).wire
    >>> f2 = FaceByPlanarWire(w2).face
    >>> builder = ShellBySewing([f1, f2])
    >>> builder.nshells
    1
    >>> len(ExploreShape.get_faces(builder.shell))
    2
    >>>
    >>> # Non-manifold case
    >>> p1 = (0., 0., 0.)
    >>> p2 = (-1., 0., 0.)
    >>> p3 = (-1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> w3 = WireByPoints([p1, p2, p3, p4], True).wire
    >>> f3 = FaceByPlanarWire(w3).face
    >>> builder = ShellBySewing([f1, f2, f3], non_manifold=True)
    >>> builder.nshells
    1
    >>> len(ExploreShape.get_faces(builder.shell))
    3
    """

    def __init__(self, faces, tol=None, cut_free_edges=False,
                 non_manifold=False):
        if tol is None:
            tol = max([ExploreShape.global_tolerance(f, 1) for f in faces])

        tool = BRepBuilderAPI_Sewing(tol, True, True, cut_free_edges,
                                     non_manifold)
        for f in faces:
            tool.Add(f)
        tool.Perform()

        shape = tool.SewedShape()
        self._shells = []
        self._shell = None
        self._nshells = 0
        if shape.ShapeType() == TopAbs_COMPOUND:
            self._shells = ExploreShape.get_shells(shape)
            self._nshells = len(self._shells)
        elif shape.ShapeType() == TopAbs_FACE:
            self._shell = CheckShape.to_shell(shape)
            self._nshells = 1
        elif shape.ShapeType() == TopAbs_SHELL:
            self._shell = TopoDS.Shell_(shape)
            self._nshells = 1
        else:
            msg = 'Invalid shape type in ShellBySewing.'
            raise RuntimeError(msg)

    @property
    def nshells(self):
        """
        :return: Number of shells.
        :rtype: int
        """
        return self._nshells

    @property
    def shell(self):
        """
        :return: The sewn shell.
        :rtype: OCCT.TopoDS.TopoDS_Shell
        """
        return self._shell

    @property
    def shells(self):
        """
        :return: The sewn shells if more than one is found.
        :rtype: list[OCCT.TopoDS.TopoDS_Shell]
        """
        return self._shells


class SolidByShell(object):
    """
    Create a solid using a shell. The shell can either be closed (finite
    solid) or open (infinite solid).

    :param OCCT.TopoDS.TopoDS_Shell shell: The shell.
    """

    def __init__(self, shell):
        self._solid = BRepBuilderAPI_MakeSolid(shell).Solid()

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._solid


# SOLID -----------------------------------------------------------------------

class SolidByPlane(object):
    """
    Create a solid box using a plane. The plane will be extruded in the
    direction of the plane's normal. The solid's width and height will be
    centered at the plane's origin.

    :param afem.geometry.entities.Plane pln: The plane.
    :param float width: Width of the box.
    :param float height: Height of the box.
    :param float depth: Depth of the box.

    Usage:

    >>> from afem.geometry import PlaneByNormal
    >>> from afem.topology import SolidByPlane
    >>> pln = PlaneByNormal().plane
    >>> box = SolidByPlane(pln, 1., 1., 1.).solid
    """

    def __init__(self, pln, width, height, depth):
        w = width / 2.
        h = height / 2.
        gp_pln = pln.handle.Pln()
        face = BRepBuilderAPI_MakeFace(gp_pln, -w, w, -h, h).Face()
        vn = pln.norm(0., 0.)
        vn.Normalize()
        vn.Scale(depth)
        shape = BRepPrimAPI_MakePrism(face, vn).Shape()
        self._solid = TopoDS.Solid_(shape)

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._solid


class SolidByDrag(object):
    """
    Create a solid by dragging a face along a vector.

    :param OCCT.TopoDS.TopoDS_Face face: The face.
    :param vector_like v: The vector to drag the shape.

    Usage:

    >>> from afem.topology import FaceByPlanarWire, SolidByDrag, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> face = FaceByPlanarWire(wire).face
    >>> builder = SolidByDrag(face, (0., 0., 1.))
    >>> solid = builder.solid
    >>> f1 = builder.first_face
    >>> f2 = builder.last_face
    """

    def __init__(self, face, v):
        v = CheckGeom.to_vector(v)
        builder = BRepPrimAPI_MakePrism(face, v)
        self._solid = TopoDS.Solid_(builder.Shape())
        self._f1 = TopoDS.Face_(builder.FirstShape())
        self._f2 = TopoDS.Face_(builder.LastShape())

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._solid

    @property
    def first_face(self):
        """
        :return: The face at the bottom of the solid.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._f1

    @property
    def last_face(self):
        """
        :return: The face at the top of the solid.
        :rtype: OCCT.TopoDS.TopoDS_Face
        """
        return self._f2


class SolidByCylinder(object):
    """
    Create a cylindrical solid.

    :param float radius: The radius.
    :param float height: The height.
    :param axis2: Not yet implemented. Solid will be constructed in xy-plane.

    :raise NotImplementedError: If an axis is provided.

    Usage:

    >>> from afem.topology import SolidByCylinder
    >>> builder = SolidByCylinder(1., 1.)
    >>> solid=builder.solid
    """

    def __init__(self, radius, height, axis2=None):
        if axis2 is None:
            builder = BRepPrimAPI_MakeCylinder(radius, height)
        else:
            raise NotImplementedError('Providing Axis2 not yet implemented.')

        self._solid = builder.Solid()

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._solid


# COMPOUND --------------------------------------------------------------------

class CompoundByShapes(object):
    """
    Create a compound from a list of shapes.

    :param list[OCCT.TopoDS.TopoDS_Shape] shapes: List of shapes.
    """

    def __init__(self, shapes):
        cp = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(cp)
        for shape in shapes:
            builder.Add(cp, shape)
        self._cp = cp

    @property
    def compound(self):
        """
        :return: The compound.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        return self._cp


# HALFSPACE -------------------------------------------------------------------

class HalfspaceByShape(object):
    """
    Create a half-space by a face or shell and a reference point.A half-space
    is an infinite solid, limited by a surface. It is built from a face or a
    shell, which bounds it, and with a reference point, which specifies the
    side of the surface where the matter of the half-space is located. A
    half-space is a tool commonly used in topological operations to cut another
    shape.

    :param shape: The face or shell.
    :type shape: OCCT.TopoDS.TopoDS_Face or OCCT.TopoDS.TopoDS_Shell
    :param point_like pnt: The reference point where the matter is located.

    :raise RuntimeError: If *shape* is not a face or shell.

    Usage:

    >>> from afem.geometry import PlaneByNormal
    >>> from afem.topology import FaceByPlane, HalfspaceByShape
    >>> pln = PlaneByNormal().plane
    >>> f = FaceByPlane(pln, -1., 1., -1., 1.).face
    >>> builder = HalfspaceByShape(f, (0., 0., 1.))
    >>> hs = builder.solid
    """

    def __init__(self, shape, pnt):
        pnt = CheckGeom.to_point(pnt)

        if shape.ShapeType() not in [TopAbs_FACE, TopAbs_SHELL]:
            msg = 'Invalid shape for creating half-space.'
            raise RuntimeError(msg)
        self._solid = BRepPrimAPI_MakeHalfSpace(shape, pnt).Solid()

    @property
    def solid(self):
        """
        :return: The half-space as a solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._solid


class HalfspaceBySurface(HalfspaceByShape):
    """
    Create a half-space using a surface.

    :param afem.geometry.entities.Surface: The surface.
    :param point_like pnt: The reference point where the matter is located.
    """

    def __init__(self, srf, pnt):
        face = FaceBySurface(srf).face
        super(HalfspaceBySurface, self).__init__(face, pnt)


# SHAPE -----------------------------------------------------------------------

class ShapeByFaces(object):
    """
    Create either a face or a shell from a list of faces.

    :param list[OCCT.TopoDS.TopoDS_Face] faces: The faces.
    :param bool sew: Option to sew the faces if more than one face is provided.
    :param float tol: Sewing tolerance. If *None* the maximum tolerance of all
        the faces is used.
    :param bool cut_free_edges: Option for cutting of free edges.
    :param bool non_manifold: Option for non-manifold processing.
    """

    def __init__(self, faces, sew=False, tol=None, cut_free_edges=False,
                 non_manifold=False):

        if len(faces) == 1:
            self._shape = faces[0]
        elif len(faces) > 1 and not sew:
            self._shape = ShellByFaces(faces).shell
        elif len(faces) > 1 and sew:
            builder = ShellBySewing(faces, tol, cut_free_edges, non_manifold)
            if builder.nshells > 1:
                self._shape = CompoundByShapes(builder.shells).compound
            elif builder.nshells == 1:
                self._shape = builder.shell
            else:
                msg = 'Failed to create any shells.'
                raise RuntimeError(msg)
        else:
            msg = 'No faces provided.'
            raise ValueError(msg)

    @property
    def shape(self):
        """
        :return: The shape. Will be a face if there was only one face in the
            list, it will be a shell if only one shell was formed, or it will
            be a compound of shells if more than one shell was found after
            sewing.
        :rtype: OCCT.TopoDS.TopoDS_Face or OCCT.TopoDS.TopoDS_Shell or
            OCCT.TopoDS.TopoDS_Compound
        """
        return self._shape

    @property
    def is_face(self):
        """
        :return: *True* if the shape is a face, *False* if not.
        :rtype: bool
        """
        return isinstance(self.shape, TopoDS_Face)

    @property
    def is_shell(self):
        """
        :return: *True* if the shape is a shell, *False* if not.
        :rtype: bool
        """
        return isinstance(self.shape, TopoDS_Shell)

    @property
    def is_compound(self):
        """
        :return: *True* if the shape is a compound, *False* if not.
        :rtype: bool
        """
        return isinstance(self.shape, TopoDS_Compound)


class ShapeByDrag(object):
    """
    Create a shape by dragging another shape along a vector.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
    :param vector_like v: The vector to drag the shape.
    """

    def __init__(self, shape, v):
        v = CheckGeom.to_vector(v)
        builder = BRepPrimAPI_MakePrism(shape, v)
        self._shape = builder.Shape()
        self._s1 = builder.FirstShape()
        self._s2 = builder.LastShape()

    @property
    def shape(self):
        """
        :return: The shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape

    @property
    def first_shape(self):
        """
        :return: The shape at the bottom.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._s1

    @property
    def last_shape(self):
        """
        :return: The shape at the top.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._s2

    @property
    def is_edge(self):
        """
        :return: *True* if shape is an edge, *False* if not.
        :rtype: bool
        """
        return self.shape.ShapeType() == TopAbs_EDGE

    @property
    def is_face(self):
        """
        :return: *True* if shape is a face, *False* if not.
        :rtype: bool
        """
        return self.shape.ShapeType() == TopAbs_FACE

    @property
    def is_shell(self):
        """
        :return: *True* if shape is a shell, *False* if not.
        :rtype: bool
        """
        return self.shape.ShapeType() == TopAbs_SHELL

    @property
    def is_solid(self):
        """
        :return: *True* if shape is a solid, *False* if not.
        :rtype: bool
        """
        return self.shape.ShapeType() == TopAbs_SOLID

    @property
    def is_compsolid(self):
        """
        :return: *True* if shape is a compsolid, *False* if not.
        :rtype: bool
        """
        return self.shape.ShapeType() == TopAbs_COMPSOLID

    @property
    def is_compound(self):
        """
        :return: *True* if shape is a compound, *False* if not.
        :rtype: bool
        """
        return self.shape.ShapeType() == TopAbs_COMPOUND


# GEOMETRY --------------------------------------------------------------------

class PointAlongShape(object):
    """
    Create a point along an edge or wire at a specified distance from the first
    parameter.

    :param shape: The shape.
    :type shape: OCCT.TopoDS.TopoDS_Edge or OCCT.TopoDS.TopoDS_Wire
    :param float ds: The distance along the curve from the given parameter.
    :param float tol: Tolerance.

    :raise TypeError: If *shape* if not a curve or wire.
    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.topology import *
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> tool = PointAlongShape(e, 5.)
    >>> tool.point
    Point(5.000, 0.000, 0.000)
    >>> tool.parameter
    5.0
    """

    def __init__(self, shape, ds, tol=1.0e-7):
        if shape.ShapeType() not in [TopAbs_EDGE, TopAbs_WIRE]:
            msg = 'Invalid shape in PointsAlongShapeByNumber.'
            raise TypeError(msg)

        shape = CheckShape.to_shape(shape)

        if shape.ShapeType() == TopAbs_EDGE:
            adp_crv = BRepAdaptor_Curve(shape)
        elif shape.ShapeType() == TopAbs_WIRE:
            adp_crv = BRepAdaptor_CompCurve(shape)
        else:
            msg = 'Could not create adaptor curve in PointAlongShape.'
            raise RuntimeError(msg)

        u0 = adp_crv.FirstParameter()
        ap = GCPnts_AbscissaPoint(tol, adp_crv, ds, u0)
        if not ap.IsDone():
            msg = "GCPnts_AbscissaPoint failed."
            raise RuntimeError(msg)

        u = ap.Parameter()
        self._u = u
        p = adp_crv.Value(u)
        self._p = Point(p.X(), p.Y(), p.Z())

    @property
    def point(self):
        """
        :return: The created point.
        :rtype: afem.geometry.entities.Point
        """
        return self._p

    @property
    def parameter(self):
        """
        :return: The parameter on the adaptor curve.
        :rtype: float
        """
        return self._u


class PointsAlongShapeByNumber(object):
    """
    Create a specified number of points along an edge or wire.

    :param shape: The shape.
    :type shape: OCCT.TopoDS.TopoDS_Edge or OCCT.TopoDS.TopoDS_Wire
    :param int n: Number of points to create (*n* > 0).
    :param float d1: An offset distance for the first point. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param OCCT.TopoDS.TopoDS_Shape shape1: A shape to define the first point.
        This shape is intersected with the edge or wire.
    :param OCCT.TopoDS.TopoDS_Shape shape2: A shape to define the last point.
        This shape is intersected with the edge or wire.
    :param float tol: Tolerance.

    :raise TypeError: If *shape* if not a curve or wire.
    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.topology import PointsAlongShapeByNumber, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> builder = PointsAlongShapeByNumber(wire, 10)
    >>> builder.npts
    10
    """

    def __init__(self, shape, n, d1=None, d2=None, shape1=None, shape2=None,
                 tol=1.0e-7):
        if shape.ShapeType() not in [TopAbs_EDGE, TopAbs_WIRE]:
            msg = 'Invalid shape in PointsAlongShapeByNumber.'
            raise TypeError(msg)

        shape1 = CheckShape.to_shape(shape1)
        shape2 = CheckShape.to_shape(shape2)

        if shape.ShapeType() == TopAbs_EDGE:
            adp_crv = BRepAdaptor_Curve(shape)
        elif shape.ShapeType() == TopAbs_WIRE:
            adp_crv = BRepAdaptor_CompCurve(shape)
        else:
            msg = 'Could not create adaptor curve in PointsAlongShapeByNumber.'
            raise RuntimeError(msg)

        u1 = adp_crv.FirstParameter()
        u2 = adp_crv.LastParameter()

        # Adjust parameters of start/end shapes are provided
        if isinstance(shape1, TopoDS_Shape):
            shape = BRepAlgoAPI_Section(shape, shape1).Shape()
            verts = ExploreShape.get_vertices(shape)
            prms = []
            for v in verts:
                gp_pnt = BRep_Tool.Pnt(v)
                proj = GeomAPI_ProjectPointOnCurve(gp_pnt, adp_crv.Curve())
                if proj.NbPoints() < 1:
                    continue
                umin = proj.LowerDistanceParameter()
                prms.append(umin)
            u1 = min(prms)

        if isinstance(shape2, TopoDS_Shape):
            shape = BRepAlgoAPI_Section(shape, shape2).Shape()
            verts = ExploreShape.get_vertices(shape)
            prms = []
            for v in verts:
                gp_pnt = BRep_Tool.Pnt(v)
                proj = GeomAPI_ProjectPointOnCurve(gp_pnt, adp_crv.Curve())
                if proj.NbPoints() < 1:
                    continue
                umin = proj.LowerDistanceParameter()
                prms.append(umin)
            u2 = max(prms)

        # Adjust u1 and u2 is d1 and d2 are provided
        if d1 is not None:
            pac = GCPnts_AbscissaPoint(tol, adp_crv, d1, u1)
            if pac.IsDone():
                u1 = pac.Parameter()
        if d2 is not None:
            pac = GCPnts_AbscissaPoint(tol, adp_crv, d2, u2)
            if pac.IsDone():
                u2 = pac.Parameter()

        # Create uniform abscissa
        ua = GCPnts_UniformAbscissa(adp_crv, int(n), u1, u2, tol)
        if not ua.IsDone():
            msg = "GCPnts_UniformAbscissa failed."
            raise RuntimeError(msg)

        self._npts = ua.NbPoints()
        self._pnts = []
        for i in range(1, self._npts + 1):
            u = ua.Parameter(i)
            p = Point()
            adp_crv.D0(u, p)
            self._pnts.append(p)

        # Point spacing
        self._ds = None
        if self._npts > 1:
            self._ds = self._pnts[0].distance(self._pnts[1])

    @property
    def npts(self):
        """
        :return: The number of points.
        :rtype: int
        """
        return self._npts

    @property
    def points(self):
        """
        :return: The points.
        :rtype: list[afem.geometry.entities.Point]
        """
        return self._pnts

    @property
    def spacing(self):
        """
        :return: The spacing between the first and second points if there
            are more than one point. Otherwise *None*.
        :rtype: float or None
        """
        return self._ds


class PointsAlongShapeByDistance(object):
    """
    Create a specified number of points along an edge or wire.

    :param shape: The shape.
    :type shape: OCCT.TopoDS.TopoDS_Edge or OCCT.TopoDS.TopoDS_Wire
    :param float maxd: The maximum allowed spacing between points. The
        actual spacing will be adjusted to not to exceed this value.
    :param float d1: An offset distance for the first point. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param OCCT.TopoDS.TopoDS_Shape shape1: A shape to define the first point.
        This shape is intersected with the edge or wire.
    :param OCCT.TopoDS.TopoDS_Shape shape2: A shape to define the last point.
        This shape is intersected with the edge or wire.
    :param int nmin: Minimum number of points to create.
    :param float tol: Tolerance.

    :raise TypeError: If *shape* if not a curve or wire.
    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.topology import PointsAlongShapeByDistance, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4], True).wire
    >>> builder = PointsAlongShapeByDistance(wire, 1.)
    >>> builder.npts
    5
    """

    def __init__(self, shape, maxd, d1=None, d2=None, shape1=None, shape2=None,
                 nmin=0, tol=1.0e-7):
        if shape.ShapeType() not in [TopAbs_EDGE, TopAbs_WIRE]:
            msg = 'Invalid shape in PointsAlongShapeByNumber.'
            raise TypeError(msg)

        shape1 = CheckShape.to_shape(shape1)
        shape2 = CheckShape.to_shape(shape2)

        if shape.ShapeType() == TopAbs_EDGE:
            adp_crv = BRepAdaptor_Curve(shape)
        elif shape.ShapeType() == TopAbs_WIRE:
            adp_crv = BRepAdaptor_CompCurve(shape)
        else:
            msg = 'Could not create adaptor curve in PointsAlongShapeByNumber.'
            raise RuntimeError(msg)

        u1 = adp_crv.FirstParameter()
        u2 = adp_crv.LastParameter()

        # Adjust parameters of start/end shapes are provided
        if isinstance(shape1, TopoDS_Shape):
            shape = BRepAlgoAPI_Section(shape, shape1).Shape()
            verts = ExploreShape.get_vertices(shape)
            prms = []
            for v in verts:
                gp_pnt = BRep_Tool.Pnt(v)
                proj = GeomAPI_ProjectPointOnCurve(gp_pnt, adp_crv.Curve())
                if proj.NbPoints() < 1:
                    continue
                umin = proj.LowerDistanceParameter()
                prms.append(umin)
            u1 = min(prms)

        if isinstance(shape2, TopoDS_Shape):
            shape = BRepAlgoAPI_Section(shape, shape2).Shape()
            verts = ExploreShape.get_vertices(shape)
            prms = []
            for v in verts:
                gp_pnt = BRep_Tool.Pnt(v)
                proj = GeomAPI_ProjectPointOnCurve(gp_pnt, adp_crv.Curve())
                if proj.NbPoints() < 1:
                    continue
                umin = proj.LowerDistanceParameter()
                prms.append(umin)
            u2 = max(prms)

        # Adjust u1 and u2 is d1 and d2 are provided
        if d1 is not None:
            pac = GCPnts_AbscissaPoint(tol, adp_crv, d1, u1)
            if pac.IsDone():
                u1 = pac.Parameter()
        if d2 is not None:
            pac = GCPnts_AbscissaPoint(tol, adp_crv, d2, u2)
            if pac.IsDone():
                u2 = pac.Parameter()

        # Determine number of points
        arc_length = GCPnts_AbscissaPoint.Length_(adp_crv, u1, u2, tol)
        n = ceil(arc_length / maxd) + 1
        if n < nmin:
            n = nmin

        # Create uniform abscissa
        ua = GCPnts_UniformAbscissa(adp_crv, int(n), u1, u2, tol)
        if not ua.IsDone():
            msg = "GCPnts_UniformAbscissa failed."
            raise RuntimeError(msg)

        self._npts = ua.NbPoints()
        self._pnts = []
        for i in range(1, self._npts + 1):
            u = ua.Parameter(i)
            p = Point()
            adp_crv.D0(u, p)
            self._pnts.append(p)

        # Point spacing
        self._ds = None
        if self._npts > 1:
            self._ds = self._pnts[0].distance(self._pnts[1])

    @property
    def npts(self):
        """
        :return: The number of points.
        :rtype: int
        """
        return self._npts

    @property
    def points(self):
        """
        :return: The points.
        :rtype: list[afem.geometry.entities.Point]
        """
        return self._pnts

    @property
    def spacing(self):
        """
        :return: The spacing between the first and second points if there
            are more than one point. Otherwise *None*.
        :rtype: float or None
        """
        return self._ds


class PlaneByEdges(object):
    """
    Create a plane by fitting it to all the edges of a shape.

    :param OCCT.TopoDS.TopoDS_Shape shape: The shape containing the edges.
    :param float tol: Edges must be within this planar tolerance. The
        tolerance is the largest value between the value provided or the
        largest tolerance of any one of the edges in the shape.

    Usage:

    >>> from afem.topology import PlaneByEdges, WireByPoints
    >>> p1 = (0., 0., 0.)
    >>> p2 = (1., 0., 0.)
    >>> p3 = (1., 1., 0.)
    >>> p4 = (0., 1., 0.)
    >>> wire = WireByPoints([p1, p2, p3, p4]).wire
    >>> builder = PlaneByEdges(wire)
    >>> builder.found
    True
    >>> pln = builder.plane
    """

    def __init__(self, shape, tol=-1.):
        builder = BRepBuilderAPI_FindPlane(shape, tol)
        self._found = builder.Found()
        self._pln = None
        if self._found:
            self._pln = Plane(builder.Plane())

    @property
    def found(self):
        """
        :return: *True* if plane was found, *False* if not.
        :rtype: bool
        """
        return self._found

    @property
    def plane(self):
        """
        :return: The plane. Returns *None* if no plane was found.
        :rtype: afem.geometry.entities.Plane
        """
        return self._pln


class PlaneByIntersectingShapes(object):
    """
    Create a plane by intersection two shapes. If no additional point is
    provided, then :class:`.PlaneByEdges`. is used. If a point is provided,
    then the edges are tessellated and the point is added to this list. Then
    the tool :class:`.PlaneByApprox` is used.

    :param OCCT.TopoDS.TopoDS_Shape shape1: The first shape.
    :param OCCT.TopoDS.TopoDS_Shape shape2: The second shape.
    :param afem.geometry.entities.Point pnt: Additional point to add to the
        edges since they might be collinear.
    :param float tol: Edges must be within this planar tolerance. The
        tolerance is the largest value between the value provided or the
        largest tolerance of any one of the edges in the shape.

    :raises ValueError: If there are less than three points after
        tessellating the edges.
    """

    def __init__(self, shape1, shape2, pnt=None, tol=-1.):
        shape = IntersectShapes(shape1, shape2).shape

        pnt = CheckGeom.to_point(pnt)

        self._found = False
        self._pln = None
        if pnt is None:
            tool = PlaneByEdges(shape, tol)
            self._found = tool.found
            self._pln = tool.plane
        elif CheckGeom.is_point(pnt):
            edges = ExploreShape.get_edges(shape)
            pnts = [pnt]
            for edge in edges:
                BRepMesh_IncrementalMesh(edge, 0.001)
                hndl_poly3d = BRep_Tool().Polygon3D(edge, edge.Location())
                if hndl_poly3d.IsNull():
                    continue
                tcol_pnts = hndl_poly3d.Nodes()
                for i in range(1, tcol_pnts.Length() + 1):
                    gp_pnt = tcol_pnts.Value(i)
                    pnt = CheckGeom.to_point(gp_pnt)
                    pnts.append(pnt)
            if len(pnts) < 3:
                msg = 'Less than three points to fit a plane.'
                raise ValueError(msg)
            if tol < 0.:
                tol = 1.0e-7
            tool = PlaneByApprox(pnts, tol)
            self._found = True
            self._pln = tool.plane
        else:
            msg = 'Invalid input.'
            raise TypeError(msg)

    @property
    def found(self):
        """
        :return: *True* if plane was found, *False* if not.
        :rtype: bool
        """
        return self._found

    @property
    def plane(self):
        """
        :return: The plane. Returns *None* if no plane was found.
        :rtype: afem.geometry.entities.Plane
        """
        return self._pln


if __name__ == "__main__":
    import doctest

    doctest.testmod()
