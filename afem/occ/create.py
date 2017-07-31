from OCC.BOPAlgo import BOPAlgo_MakerVolume
from OCC.BRep import BRep_Builder, BRep_Tool
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.BRepAlgo import brepalgo_ConcatenateWireC0
from OCC.BRepBuilderAPI import BRepBuilderAPI_Copy, BRepBuilderAPI_MakeEdge, \
    BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCC.BRepMesh import BRepMesh_IncrementalMesh
from OCC.BRepOffset import BRepOffset_Pipe, BRepOffset_RectoVerso, \
    BRepOffset_Skin
from OCC.BRepOffsetAPI import BRepOffsetAPI_MakeOffset, \
    BRepOffsetAPI_MakeOffsetShape, BRepOffsetAPI_MakePipeShell, \
    BRepOffsetAPI_NormalProjection
from OCC.BRepPrimAPI import BRepPrimAPI_MakeHalfSpace, BRepPrimAPI_MakePrism
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_UniformAbscissa
from OCC.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.GeomAbs import GeomAbs_Arc, GeomAbs_Intersection, GeomAbs_Tangent
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.ShapeAnalysis import ShapeAnalysis_FreeBounds_ConnectEdgesToWires
from OCC.TopAbs import TopAbs_EDGE, \
    TopAbs_REVERSED, TopAbs_WIRE
from OCC.TopTools import Handle_TopTools_HSequenceOfShape, \
    TopTools_HSequenceOfShape
from OCC.TopoDS import TopoDS_Compound, TopoDS_Face, TopoDS_Shape, \
    TopoDS_Shell, TopoDS_Wire, topods_Wire

from .check import CheckShape
from .explore import ExploreShape
from ..geometry import BBox, CheckGeom, CreateGeom


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


class ShapeBuilder(object):
    """
    Base class for building shapes.

    :var bool success:
    :var TopoDS_Shape shape: Created shape.
    :var list shapes: List of created shapes if applicable.
    """

    def __init__(self):
        self._success = False
        self._performed = False
        self._shape = None
        self._shapes = []

    @property
    def success(self):
        if not self._performed:
            self.build()
        return self._success

    @property
    def shape(self):
        if not self._performed:
            self.build()
        return self._shape

    @property
    def shapes(self):
        if not self._performed:
            self.build()
        return self._shapes

    def build(self):
        """
        Build shape. Should be overridden.
        """
        pass


class CompoundByShapes(ShapeBuilder):
    """
    Create a compound from a list of shapes.
    """

    def __init__(self, shapes):
        super(CompoundByShapes, self).__init__()
        self._shapes_in = []
        for shape in shapes:
            shape = CheckShape.to_shape(shape)
            if shape:
                self._shapes_in.append(shape)

    def build(self):
        self._performed = True
        cp = TopoDS_Compound()
        builder = BRep_Builder()
        builder.MakeCompound(cp)
        for shape in self._shapes_in:
            builder.Add(cp, shape)
        self._shape = cp
        self._success = True

    @property
    def compound(self):
        return self.shape


class BoxByPlane(object):
    pass


class WireByConnectedEdges(object):
    pass


class ShapeBySweep(object):
    pass


class EdgeBySweep(object):
    pass


class FaceBySweep(object):
    pass


class ShellBySweep(object):
    pass


class SolidBySweep(object):
    pass


class CompSolidBySweep(object):
    pass


class WireByConcat(object):
    pass


class PlaneBySection(object):
    pass


class PointsAlongShape(object):
    pass


class PointsAlongEdge(object):
    pass


class PointsAlongWire(object):
    pass


class FaceByPlane(object):
    pass


class WireByOffset(object):
    pass


class ShellByPipe(object):
    pass


class HalfspaceByShape(object):
    pass


class WiresFromShape(object):
    pass


class ShapeByOffset(object):
    pass


class SolidsFromShapes(object):
    pass
