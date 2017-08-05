from OCC.BOPAlgo import BOPAlgo_MakerVolume
from OCC.BRep import BRep_Builder, BRep_Tool
from OCC.BRepAdaptor import BRepAdaptor_CompCurve, BRepAdaptor_Curve
from OCC.BRepAlgo import brepalgo_ConcatenateWireC0
from OCC.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.BRepBuilderAPI import BRepBuilderAPI_Copy, BRepBuilderAPI_FindPlane, \
    BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace, \
    BRepBuilderAPI_MakePolygon, BRepBuilderAPI_MakeShell, \
    BRepBuilderAPI_MakeSolid, BRepBuilderAPI_MakeVertex, \
    BRepBuilderAPI_MakeWire, BRepBuilderAPI_Sewing
from OCC.BRepMesh import BRepMesh_IncrementalMesh
from OCC.BRepOffset import BRepOffset_Pipe, BRepOffset_RectoVerso, \
    BRepOffset_Skin
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset, \
    BRepOffsetAPI_MakeOffsetShape, BRepOffsetAPI_MakePipeShell, \
    BRepOffsetAPI_NormalProjection
from OCC.BRepPrimAPI import BRepPrimAPI_MakeHalfSpace, BRepPrimAPI_MakePrism
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
from OCC.GEOMAlgo import GEOMAlgo_Splitter
from OCC.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.GeomAbs import GeomAbs_Arc, GeomAbs_Intersection, GeomAbs_Tangent
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBounds_ConnectEdgesToWires
from OCC.ShapeBuild import ShapeBuild_ReShape
from OCC.TopAbs import TopAbs_COMPOUND, TopAbs_EDGE, TopAbs_FACE, \
    TopAbs_REVERSED, TopAbs_SHELL, TopAbs_WIRE
from OCC.TopTools import Handle_TopTools_HSequenceOfShape, \
    TopTools_HSequenceOfShape
from OCC.TopoDS import TopoDS_Compound, TopoDS_Face, TopoDS_Shape, \
    TopoDS_Shell, TopoDS_Wire, topods_Edge, topods_Face, topods_Shell, \
    topods_Solid, topods_Vertex, topods_Wire
from numpy import ceil

from afem.geometry import BBox, CheckGeom, CreateGeom
from afem.geometry.entities import Plane, Point
from afem.topology.check import CheckShape
from afem.topology.explore import ExploreShape
from afem.topology.tools import ShapeTools

__all__ = ["VertexByPoint", "EdgeByPoints", "EdgeByVertices", "EdgeByCurve",
           "EdgeByDrag", "WireByEdges", "WiresByConnectedEdges",
           "WireByConcat",
           "WireByPlanarOffset", "WiresByShape", "WireByPoints", "WireBySplit",
           "FaceBySurface",
           "FaceByPlane", "FaceByPlanarWire", "FaceByDrag", "ShellBySurface",
           "ShellByFaces", "ShellBySewing",
           "SolidByShell", "ShellByDrag", "SolidByPlane", "SolidByDrag",
           "CompoundByShapes", "HalfspaceByShape", "PointsAlongShapeByNumber",
           "PointsAlongShapeByDistance", "PlaneByEdges"]

_occ_join = {'a': GeomAbs_Arc,
             'arc': GeomAbs_Arc,
             't': GeomAbs_Tangent,
             'tangent': GeomAbs_Tangent,
             'i': GeomAbs_Intersection,
             'intersection': GeomAbs_Intersection}


class CreateShape(object):
    """
    Shape creator.
    """

    @staticmethod
    def shape_by_sweep(shape, vec):
        """
        Make a linear swept topology (i.e., prism).

        :param shape:
        :param vec:

        :return:
        """
        shape = CheckShape.to_shape(shape)
        vec = CheckGeom.to_vector(vec)
        return BRepPrimAPI_MakePrism(shape, vec).Shape()

    @staticmethod
    def shape_by_offset(shape, distance, tol, mode=BRepOffset_Skin,
                        intersection=False, self_intersect=False,
                        join=GeomAbs_Arc):
        """
        Offset a shape.

        :param shape:
        :param distance:
        :param tol:
        :param mode:
        :param intersection:
        :param self_intersect:
        :param join:

        :return:
        """
        shape = CheckShape.to_shape(shape)
        if not shape:
            return None

        if join < GeomAbs_Arc:
            join = GeomAbs_Arc
        if join > GeomAbs_Intersection:
            join = GeomAbs_Intersection
        if join not in [GeomAbs_Arc, GeomAbs_Tangent, GeomAbs_Intersection]:
            return None

        if mode < BRepOffset_Skin:
            mode = BRepOffset_Skin
        if mode > BRepOffset_RectoVerso:
            mode = BRepOffset_RectoVerso
        if mode not in [BRepOffset_Skin, BRepOffset_Pipe,
                        BRepOffset_RectoVerso]:
            return None

        offset = BRepOffsetAPI_MakeOffsetShape(shape, distance, tol, mode,
                                               intersection, self_intersect,
                                               join)
        if not offset.IsDone():
            return None

        return CheckShape.to_shape(offset.Shape())

    @staticmethod
    def shape_by_pipe(spine, profiles, spine_support,
                      with_contact=False, with_correction=False):
        """
        Create a shell by sweeping profile(s) along a spine with a normal set
        by a shape.

        :param spine:
        :param profiles:
        :param spine_support:
        :param with_contact:
        :param with_correction:

        :return:
        """
        spine = CheckShape.to_wire(spine)
        spine_support = CheckShape.to_shape(spine_support)
        if not spine or not spine_support:
            return None

        builder = BRepOffsetAPI_MakePipeShell(spine)
        builder.SetMode(spine_support)
        for profile in profiles:
            profile = CheckShape.to_shape(profile)
            if not profile:
                continue
            builder.Add(profile, with_contact, with_correction)
        if not builder.IsReady():
            return None
        builder.Build()
        if not builder.IsDone():
            return None
        return CheckShape.to_shape(builder.Shape())

    @staticmethod
    def volumes_by_shapes(shapes, intersect=False):
        """
        Make volume(s) from a list of shapes.

        :param list shapes: List of shapes to create volume(s).
        :param bool intersect: Option to first intersect the shapes.

        :return: Volume(s) created from list of shapes. If a single volume
            is created, a solid will be returned. If more than one
            volume was created, a compound will be returned containing all
            the solids. If nothing is created *None* will be returned.
        :rtype: TopoDS_Shape or None
        """
        bop = BOPAlgo_MakerVolume()
        for shape in shapes:
            shape = CheckShape.to_shape(shape)
            if not shape:
                continue
            bop.AddArgument(shape)
        bop.SetIntersect(intersect)
        bop.Perform()

        if bop.ErrorStatus() != 0:
            return None

        return bop.Shape()

    @staticmethod
    def wire_by_concat(wire):
        """
        Create a wire by joining all the edges into a single wire. The
        wire may have C0 continuity.

        :param wire:
        :return:
        """
        wire = CheckShape.to_wire(wire)
        if not wire:
            return None

        return brepalgo_ConcatenateWireC0(wire)

    @staticmethod
    def wires_by_edges(edges, tol=None, shared=False):
        """
        Build wires from a list of unsorted edges.

        :param edges:
        :param tol:
        :param shared:

        :return:
        """
        hedges = TopTools_HSequenceOfShape()
        for e in edges:
            hedges.Append(e)

        if tol is None:
            tol = max([ExploreShape.get_tolerance(e, 1) for e in edges])

        # noinspection PyArgumentList
        hwires = Handle_TopTools_HSequenceOfShape()
        ShapeAnalysis_FreeBounds_ConnectEdgesToWires(hedges.GetHandle(),
                                                     tol, shared, hwires)

        wires_obj = hwires.GetObject()
        wires = []
        for i in range(1, wires_obj.Length() + 1):
            w = topods_Wire(wires_obj.Value(i))
            wires.append(w)

        return wires

    @staticmethod
    def wires_from_shape(shape):
        """
        Use edges from the shape to create connected wires.

        :param shape:

        :return:
        """
        edges = ExploreShape.get_edges(shape)
        if not edges:
            return []

        return CreateShape.wires_by_edges(edges)

    @staticmethod
    def wire_by_offset(spine, distance, altitude=0., join=GeomAbs_Arc,
                       is_open=False):
        """
        Offset wire in a planar face.

        :param spine:
        :param distance:
        :param altitude:
        :param join:
        :param is_open:

        :return:
        """
        spine = CheckShape.to_shape(spine)
        if not isinstance(spine, (TopoDS_Wire, TopoDS_Face)):
            return None

        if join < GeomAbs_Arc:
            join = GeomAbs_Arc
        if join > GeomAbs_Intersection:
            join = GeomAbs_Intersection
        if join not in [GeomAbs_Arc, GeomAbs_Tangent, GeomAbs_Intersection]:
            return None

        offset = BRepOffsetAPI_MakeOffset(spine, join, is_open)
        offset.Perform(distance, altitude)
        if not offset.IsDone():
            return None

        return CheckShape.to_shape(offset.Shape())

    @staticmethod
    def wires_by_projection(basis_shape, to_project):
        """
        Perform the normal projection of a shape.

        :param basis_shape:
        :param to_project:
        :return:
        """
        basis_shape = CheckShape.to_shape(basis_shape)
        to_project = CheckShape.to_shape(to_project)
        if not basis_shape or not to_project:
            return []

        proj = BRepOffsetAPI_NormalProjection(basis_shape)
        proj.Add(to_project)
        proj.Build()
        if not proj.IsDone():
            return []

        shape = proj.Shape()
        return ExploreShape.get_wires(shape)

    @staticmethod
    def face_from_plane(pln, umin, umax, vmin, vmax):
        """
        Create a finite face from a plane.

        :param pln:
        :param umin:
        :param umax:
        :param vmin:
        :param vmax:

        :return:
        """
        if not CheckGeom.is_plane(pln):
            return None
        builder = BRepBuilderAPI_MakeFace(pln.Pln(), umin, umax, vmin, vmax)
        if not builder.IsDone():
            return None
        return builder.Face()

    @staticmethod
    def sew_faces(faces, tol=None):
        """
        Sew faces to make a shell or compound of shells.

        :param faces:
        :param tol:

        :return:
        """
        if tol is None:
            tol = max([ExploreShape.get_tolerance(f, 1) for f in faces])

        shell = TopoDS_Shell()
        builder = BRep_Builder()
        builder.MakeShell(shell)
        for f in faces:
            builder.Add(shell, f)

        if len(faces) == 1:
            return shell

        sew = BRepBuilderAPI_Sewing(tol)
        sew.Load(shell)
        sew.Perform()
        return sew.SewedShape()

    @staticmethod
    def box_from_plane(pln, width, height, depth):
        """
        Create a solid box from a plane.

        :param pln:
        :param width:
        :param height:
        :param depth:

        :return:
        """
        if not CheckGeom.is_plane(pln):
            return None

        # Make finite face from plane.
        w = width / 2.
        h = height / 2.
        gp_pln = pln.Pln()
        face = BRepBuilderAPI_MakeFace(gp_pln, -w, w, -h, h).Face()

        # Get normal vector at center of plane and use depth to extrude.
        vn = pln.norm(0., 0.)
        vn.Normalize()
        vn.Scale(depth)
        return BRepPrimAPI_MakePrism(face, vn).Shape()

    @staticmethod
    def halfspace_by_point(shape, pref):
        """
        Create a half-space by a shape and a reference point.

        :param shape:
        :param pref:

        :return:
        """
        shape = CheckShape.to_shape(shape)
        if not shape:
            return None

        pref = CheckGeom.to_point(pref)
        if not CheckGeom.is_point(pref):
            return None

        return BRepPrimAPI_MakeHalfSpace(shape, pref).Solid()

    @staticmethod
    def compound_by_shapes(shapes):
        """
        Create a compound from a list of shapes.

        :param shapes:

        :return:
        """
        compound = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(compound)
        for shape in shapes:
            shape = CheckShape.to_shape(shape)
            if shape:
                builder.Add(compound, shape)
        return compound

    @staticmethod
    def plane_from_section(shape1, shape2, pref):
        """
        Create plane from intersection of two shapes.
        """
        # Intersect the two shapes.
        edges = ShapeTools.bsection(shape1, shape2, 'edge')
        if not edges:
            return None

        # Tessellate the edges.
        for e in edges:
            ShapeTools.incremental_mesh(e, 1., True)

        # Gather points to fit a plane.
        pnts = [pref]
        for e in edges:
            hpoly = BRep_Tool().Polygon3D(e, e.Location())
            tcol_pnts = hpoly.GetObject().Nodes()
            for i in range(1, tcol_pnts.Length() + 1):
                gp_pnt = tcol_pnts.Value(i)
                pnt = CheckGeom.to_point(gp_pnt)
                pnts.append(pnt)
        if len(pnts) < 3:
            return None

        return CreateGeom.plane_by_fit(pnts)

    @staticmethod
    def points_along_edge(edge, dx, tol=1.0e-7):
        """
        Create points along an edge.

        :param edge:
        :param dx:
        :param tol:

        :return:
        """
        edge = CheckShape.to_edge(edge)
        if not edge:
            return []

        # Calculate number of points.
        adp_crv = BRepAdaptor_Curve(edge)
        arc_len = GCPnts_AbscissaPoint.Length(adp_crv,
                                              adp_crv.FirstParameter(),
                                              adp_crv.LastParameter(), tol)
        nb_pts = int(arc_len / dx) + 1

        pac = GCPnts_UniformAbscissa(adp_crv, nb_pts, tol)
        if not pac.IsDone():
            return []

        pnts = []
        for i in range(1, pac.NbPoints() + 1):
            u = pac.Parameter(i)
            gp_pnt = adp_crv.Value(u)
            pnts.append(CreateGeom.point_by_xyz(gp_pnt.X(), gp_pnt.Y(),
                                                gp_pnt.Z()))

        if edge.Orientation() == TopAbs_REVERSED:
            pnts.reverse()

        return pnts

    @staticmethod
    def points_along_shape(shape, maxd=None, npts=None, u1=None, u2=None,
                           s1=None, s2=None, shape1=None, shape2=None,
                           tol=1.0e-7):
        """
        Create points along an edge, wire, or curve.

        :param shape:
        :param maxd:
        :param npts:
        :param u1:
        :param u2:
        :param s1:
        :param s2:
        :param shape1:
        :param shape2:
        :param tol:

        :return:
        """
        # Get adaptor curve.
        if CheckGeom.is_curve_like(shape):
            adp_crv = GeomAdaptor_Curve(shape.GetHandle())
            edge = BRepBuilderAPI_MakeEdge(shape.GetHandle()).Edge()
        elif isinstance(shape, TopoDS_Shape):
            if shape.ShapeType() == TopAbs_EDGE:
                edge = shape
                adp_crv = BRepAdaptor_Curve(edge)
            elif shape.ShapeType() == TopAbs_WIRE:
                edge = brepalgo_ConcatenateWireC0(shape)
                adp_crv = BRepAdaptor_Curve(edge)
            else:
                return []
        else:
            return []

        # Check parameters.
        if u1 is None:
            u1 = adp_crv.FirstParameter()
        if u2 is None:
            u2 = adp_crv.LastParameter()
        if u1 > u2:
            u1, u2 = u2, u1

        # Adjust parameter if shapes are provided.
        if isinstance(shape1, TopoDS_Shape):
            verts = ShapeTools.bsection(edge, shape1, 'vertex')
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
            verts = ShapeTools.bsection(edge, shape2, 'vertex')
            prms = []
            for v in verts:
                gp_pnt = BRep_Tool.Pnt(v)
                proj = GeomAPI_ProjectPointOnCurve(gp_pnt, adp_crv.Curve())
                if proj.NbPoints() < 1:
                    continue
                umin = proj.LowerDistanceParameter()
                prms.append(umin)
            u2 = max(prms)

        # Adjust u1 and u2 is s1 and s2 are provided.
        if s1 is not None:
            pac = GCPnts_AbscissaPoint(tol, adp_crv, s1, u1)
            if pac.IsDone():
                u1 = pac.Parameter()
        if s2 is not None:
            pac = GCPnts_AbscissaPoint(tol, adp_crv, s2, u2)
            if pac.IsDone():
                u2 = pac.Parameter()

        # Adjust step size if necessary.
        if maxd is not None:
            arc_len = GCPnts_AbscissaPoint.Length(adp_crv, u1, u2, tol)
            nb_pts = int(arc_len / maxd) + 1
        else:
            nb_pts = int(npts)

        # Minimum number of points if maxd and npts are provided.
        if maxd is not None and npts is not None:
            if nb_pts < npts:
                nb_pts = int(npts)

        if nb_pts < 1:
            return []

        # OCC uniform abscissa.
        occ_pnts = GCPnts_UniformAbscissa(adp_crv, nb_pts, u1, u2, tol)
        if not occ_pnts.IsDone():
            return []

        # Gather results.
        npts = occ_pnts.NbPoints()
        points = []
        for i in range(1, npts + 1):
            u = occ_pnts.Parameter(i)
            gp_pnt = adp_crv.Value(u)
            pnt = CheckGeom.to_point(gp_pnt)
            points.append(pnt)

        return points

    @staticmethod
    def copy_shape(shape, copy_geom=True):
        """
        Copy a shape.

        :param shape:
        :param copy_geom:

        :return:
        """
        copy = BRepBuilderAPI_Copy(shape, copy_geom)
        if not copy.IsDone():
            return None
        return copy.Shape()

    @staticmethod
    def bounding_box(shape):
        """
        Create a bounding box from the shape.

        :param shape:

        :return:
        """
        shape = ExploreShape.to_shape(shape)
        if not shape:
            return None

        if shape.IsNull():
            return None

        bbox = BBox()
        bbox.add_shape(shape)
        return bbox

    @staticmethod
    def incremental_mesh(shape, linear, is_relative=False, angular=0.5):
        """
        Builds the mesh of a shape.

        :param shape:
        :param linear:
        :param is_relative:
        :param angular:

        :return:
        """
        BRepMesh_IncrementalMesh(shape, linear, is_relative, angular, False)


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
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v


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
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._e

    @property
    def vertex1(self):
        """
        :return: The first vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v1

    @property
    def vertex2(self):
        """
        :return: The second vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v2


class EdgeByVertices(object):
    """
    Create an edge between two vertices.

    :param OCC.TopoDS.TopoDS_Vertex v1: The first vertex.
    :param OCC.TopoDS.TopoDS_Vertex v2: The second vertex.

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
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._e


class EdgeByCurve(object):
    """
    Create an edge using a curve.

    :param curve_like crv: The curve.

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
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._e

    @property
    def vertex1(self):
        """
        :return: The first vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v1

    @property
    def vertex2(self):
        """
        :return: The second vertex.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v2


class EdgeByDrag(object):
    """
    Create an edge by dragging a vertex along a vector.

    :param OCC.TopoDS.TopoDS_Vertex vertex: The vertex.
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
        self._e = topods_Edge(builder.Shape())
        self._v1 = topods_Vertex(builder.FirstShape())
        self._v2 = topods_Vertex(builder.LastShape())

    @property
    def edge(self):
        """
        :return: The edge.
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._e

    @property
    def first_vertex(self):
        """
        :return: The vertex at the bottom of the edge.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v1

    @property
    def last_vertex(self):
        """
        :return: The vertex at the top of the edge.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._v2


class WireByEdges(object):
    """
    Create a wire using up to four edges.

    :param OCC.TopoDS.TopoDS_Edge e1: The first edge.
    :param OCC.TopoDS.TopoDS_Edge e2: The second edge.
    :param OCC.TopoDS.TopoDS_Edge e3: The third edge.
    :param OCC.TopoDS.TopoDS_Edge e4: The fourth edge.

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
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._w

    @property
    def last_edge(self):
        """
        :return: The last edge added to the wire.
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._last_e

    @property
    def last_vertex(self):
        """
        :return: The last vertex added to the wire.
        :rtype: OCC.TopoDS.TopoDS_Vertex
        """
        return self._last_v


class WiresByConnectedEdges(object):
    """
    Create wires from a list of unsorted edges.

    :param list[OCC.TopoDS.TopoDS_Edge] edges: List of edges.
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
            tol = max([ExploreShape.get_tolerance(e, 1) for e in edges])

        hwires = Handle_TopTools_HSequenceOfShape()
        ShapeAnalysis_FreeBounds_ConnectEdgesToWires(hedges.GetHandle(),
                                                     tol, shared, hwires)

        wires_obj = hwires.GetObject()
        wires = []
        for i in range(1, wires_obj.Length() + 1):
            w = topods_Wire(wires_obj.Value(i))
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
        :rtype: list[OCC.TopoDS.TopoDS_Wire]
        """
        return self._wires


class WireByConcat(object):
    """
    Create a wire by concatenating all the edges into a single edge. The
    wire may have C0 continuity.

    :param OCC.TopoDS.TopoDS_Wire wire: The wire.

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.topology import EdgeByCurve, WireByEdges, WireByConcat
    >>> c1 = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> c2 = NurbsCurveByPoints([(10., 0., 0.), (11., 1., 0.)]).curve
    >>> e1 =  EdgeByCurve(c1).edge
    >>> e2 =  EdgeByCurve(c2).edge
    >>> w = WireByEdges(e1, e2).wire
    >>> wnew = WireByConcat(w).wire
    """

    def __init__(self, wire):
        self._wire = brepalgo_ConcatenateWireC0(wire)

    @property
    def wire(self):
        """
        :return: The wire.
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._wire


class WireByPlanarOffset(object):
    """
    Create a wire by offsetting a planar wire or face.

    :param spine: The wire to offset. If a face is provided the outer wire
        will be used.
    :type spine: OCC.TopoDS.TopoDS_Wire or OCC.TopoDS.TopoDS_Face
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
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._w


class WiresByShape(WiresByConnectedEdges):
    """
    Create wires by connecting all the edges of a shape. This method gathers
    all the unique edges of a shape and then uses
    :class:`.WiresByConnectedEdges`.

    :param OCC.TopoDS.TopoDS_Shape: The shape.

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
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._w


class WireBySplit(object):
    """
    Split a wire with a shape.

    :param OCC.TopoDS.TopoDS_Wire wire: The wire.
    :param OCC.TopoDS.TopoDS_Shape splitter: The splitter shape.
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
            self._wire = topods_Wire(new_wire)

    @property
    def wire(self):
        """
        :return: The split wire.
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._wire


class FaceBySurface(object):
    """
    Create a face from a surface.

    :param surface_like srf: The surface.
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
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._f


class FaceByPlane(object):
    """
    Create a finite face from a plane.

    :param afem.geometry.entities.Plane: The plane.
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
        builder = BRepBuilderAPI_MakeFace(pln.Pln(), umin, umax, vmin, vmax)
        self._f = builder.Face()

    @property
    def face(self):
        """
        :return: The face.
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._f


class FaceByPlanarWire(object):
    """
    Create a face from a planar wire.

    :param OCC.TopoDS.TopoDS_Wire wire: The wire.

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
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._f


class FaceByDrag(object):
    """
    Create a face by dragging an edge along a vector.

    :param OCC.TopoDS.TopoDS_Edge edge: The edge.
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
        self._f = topods_Face(builder.Shape())
        self._e1 = topods_Edge(builder.FirstShape())
        self._e2 = topods_Edge(builder.LastShape())

    @property
    def face(self):
        """
        :return: The face.
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._f

    @property
    def first_edge(self):
        """
        :return: The edge at the bottom of the face.
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._e1

    @property
    def last_edge(self):
        """
        :return: The edge at the top of the face.
        :rtype: OCC.TopoDS.TopoDS_Edge
        """
        return self._e2


class ShellBySurface(object):
    """
    Create a shell from a surface.

    :param surface_like srf: The surface.

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
        :rtype: OCC.TopoDS.TopoDS_Shell
        """
        return self._shell


class ShellByFaces(object):
    """
    Create a shell from connected faces. This method initializes a shell an
    then simply adds the faces to it. The faces should already have shared
    edges. This is not checked.

    :param list[OCC.TopoDS.TopoDS_Face] faces: List of faces.

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
        :rtype: OCC.TopoDS.TopoDS_Shell
        """
        return self._shell


class ShellByDrag(object):
    """
    Create a shell by dragging a wire along a vector.

    :param OCC.TopoDS.TopoDS_Wire wire: The wire.
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
        self._shell = topods_Shell(builder.Shape())
        self._w1 = topods_Wire(builder.FirstShape())
        self._w2 = topods_Wire(builder.LastShape())

    @property
    def shell(self):
        """
        :return: The shell.
        :rtype: OCC.TopoDS.TopoDS_Shell
        """
        return self._shell

    @property
    def first_wire(self):
        """
        :return: The wire at the bottom of the shell.
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._w1

    @property
    def last_wire(self):
        """
        :return: The wire at the top of the shell.
        :rtype: OCC.TopoDS.TopoDS_Wire
        """
        return self._w2


class ShellBySewing(object):
    """
    Create a shell by sewing faces.

    :param list[OCC.TopoDS.TopoDS_Face] faces: The faces.
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
            tol = max([ExploreShape.get_tolerance(f, 1) for f in faces])

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
            self._shell = topods_Shell(shape)
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
        :rtype: OCC.TopoDS.TopoDS_Shell
        """
        return self._shell

    @property
    def shells(self):
        """
        :return: The sewn shells if more than one is found.
        :rtype: list[OCC.TopoDS.TopoDS_Shell]
        """
        return self._shells


class SolidByShell(object):
    """
    Create a solid using a shell. The shell can either be closed (finite
    solid) or open (infinite solid).

    :param OCC.TopoDS.TopoDS_Shell shell: The shell.
    """

    def __init__(self, shell):
        self._solid = BRepBuilderAPI_MakeSolid(shell).Solid()

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCC.TopoDS.TopoDS_Solid
        """
        return self._solid


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
        gp_pln = pln.Pln()
        face = BRepBuilderAPI_MakeFace(gp_pln, -w, w, -h, h).Face()
        vn = pln.norm(0., 0.)
        vn.Normalize()
        vn.Scale(depth)
        shape = BRepPrimAPI_MakePrism(face, vn).Shape()
        self._solid = topods_Solid(shape)

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCC.TopoDS.TopoDS_Solid
        """
        return self._solid


class SolidByDrag(object):
    """
    Create a solid by dragging a face along a vector.

    :param OCC.TopoDS.TopoDS_Face face: The face.
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
        self._solid = topods_Solid(builder.Shape())
        self._f1 = topods_Face(builder.FirstShape())
        self._f2 = topods_Face(builder.LastShape())

    @property
    def solid(self):
        """
        :return: The solid.
        :rtype: OCC.TopoDS.TopoDS_Solid
        """
        return self._solid

    @property
    def first_face(self):
        """
        :return: The face at the bottom of the solid.
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._f1

    @property
    def last_face(self):
        """
        :return: The face at the top of the solid.
        :rtype: OCC.TopoDS.TopoDS_Face
        """
        return self._f2


class CompoundByShapes(object):
    """
    Create a compound from a list of shapes.

    :param list[OCC.TopoDS.TopoDS_Shape] shapes: List of shapes.
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
        :rtype: OCC.TopoDS.TopoDS_Compound
        """
        return self._cp


class HalfspaceByShape(object):
    """
    Create a half-space by a face or shell and a reference point.A half-space
    is an infinite solid, limited by a surface. It is built from a face or a
    shell, which bounds it, and with a reference point, which specifies the
    side of the surface where the matter of the half-space is located. A
    half-space is a tool commonly used in topological operations to cut another
    shape.

    :param shape: The face or shell.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
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
        :rtype: OCC.TopoDS.TopoDS_Solid
        """
        return self._solid


class PointsAlongShapeByNumber(object):
    """
    Create a specified number of points along an edge or wire.

    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param int n: Number of points to create (*n* > 0).
    :param float d1: An offset distance for the first point. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param OCC.TopoDS.TopoDS_Shape shape1: A shape to define the first point.
        This shape is intersected with the edge or wire.
    :param OCC.TopoDS.TopoDS_Shape shape2: A shape to define the last point.
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
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param float maxd: The maximum allowed spacing between points. The
        actual spacing will be adjusted to not to exceed this value.
    :param float d1: An offset distance for the first point. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last point. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param OCC.TopoDS.TopoDS_Shape shape1: A shape to define the first point.
        This shape is intersected with the edge or wire.
    :param OCC.TopoDS.TopoDS_Shape shape2: A shape to define the last point.
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
        arc_length = GCPnts_AbscissaPoint.Length(adp_crv, u1, u2, tol)
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

    :param OCC.TopoDS.TopoDS_Shape shape: The shape containing the edges.
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
            gp_pln = builder.Plane().GetObject().Pln()
            self._pln = Plane(gp_pln)

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
