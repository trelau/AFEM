from OCC.SMESH import SMESH_Mesh
from OCC.TopoDS import TopoDS_Shape

from afem.fem.elements import Elm2D
from afem.fem.meshes import MeshData
from afem.geometry import CreateGeom
from afem.geometry.check import CheckGeom
from afem.geometry.create import PointFromParameter
from afem.geometry.project import ProjectPointToCurve, ProjectPointToSurface
from afem.graphics.viewer import ViewableItem
from afem.structure.methods.cut_parts import cut_part, \
    cut_wing_part_with_circle
from afem.structure.methods.explore_parts import get_shared_edges, \
    get_shared_nodes
from afem.structure.methods.fuse_parts import fuse_surface_part
from afem.structure.methods.merge_parts import merge_surface_part
from afem.structure.methods.modify_parts import discard_faces_by_distance, \
    discard_faces_by_solid, discard_wing_part_faces, unify_surface_part
from afem.structure.methods.reshape_parts import reshape_parts
from afem.structure.methods.sew_parts import sew_surface_parts
from afem.structure.methods.split_parts import split_part
from afem.topology import ShapeTools
from afem.topology.check import CheckShape
from afem.topology.create import PointsAlongShapeByDistance, \
    PointsAlongShapeByNumber
from afem.topology.explore import ExploreShape
from afem.topology.modify import FixShape
from afem.topology.props import LinearProps, SurfaceProps

__all__ = ["Part", "CurvePart", "Beam", "SurfacePart", "WingPart", "Spar",
           "Rib", "FuselagePart", "Bulkhead", "Floor", "Frame", "Skin",
           "Stiffener1D", "Stiffener2D", "Stringer"]


class Part(TopoDS_Shape, ViewableItem):
    """
    Base class for all parts.

    :param str label: The label.
    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """
    _indx = 1
    _all = {}

    def __init__(self, label, shape, cref=None, sref=None):
        # Initialize
        super(Part, self).__init__()
        ViewableItem.__init__(self)
        self._cref, self._sref = None, None
        self._metadata = {}
        self._subparts = {}

        # Set ID
        self._id = Part._indx
        Part._all[self._id] = self
        Part._indx += 1

        # Set attributes
        self._label = label
        self.set_shape(shape)
        self.set_cref(cref)
        self.set_sref(sref)

        print('Creating part: ', label)

    @property
    def label(self):
        """
        :return: The part label.
        :rtype: str
        """
        return self._label

    @property
    def id(self):
        """
        :return: The unique part ID.
        :rtype: int
        """
        return self._id

    @property
    def shape(self):
        """
        :return: The part shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self

    @property
    def is_null(self):
        """
        :return: *True* if part shape is null, *False* if not.
        :rtype: bool
        """
        return self.IsNull()

    @property
    def tol(self):
        """
        :return: The average tolerance of the part shape.
        :rtype: float
        """
        return ExploreShape.get_tolerance(self, 0)

    @property
    def max_tol(self):
        """
        :return: The maximum tolerance of the part shape.
        :rtype: float
        """
        return ExploreShape.get_tolerance(self, 1)

    @property
    def min_tol(self):
        """
        :return: The minimum tolerance of the part shape.
        :rtype: float
        """
        return ExploreShape.get_tolerance(self, 2)

    @property
    def metadata(self):
        """
        :return: The part metadata.
        :rtype: dict
        """
        return self._metadata

    @property
    def subparts(self):
        """
        :return: List of sub-parts associated to this part.
        :rtype: list[afem.structure.entities.Part]
        """
        return self._subparts.values()

    @property
    def cref(self):
        """
        :return: The part reference curve.
        :rtype: afem.geometry.entities.Curve
        """
        return self._cref

    @property
    def sref(self):
        """
        :return: The part reference surface.
        :rtype: afem.geometry.entities.Surface
        """
        return self._sref

    @property
    def has_cref(self):
        """
        :return: *True* if part has a reference curve, *False* if not.
        :rtype: bool
        """
        return CheckGeom.is_curve_like(self._cref)

    @property
    def has_sref(self):
        """
        :return: *True* if part has a reference surface, *False* if not.
        :rtype: bool
        """
        return CheckGeom.is_surface_like(self._sref)

    @property
    def p1(self):
        """
        :return: The first point of the part reference curve.
        :rtype: afem.geometry.entities.Point

        :raises ValueError: If the part has no reference curve.
        """
        if not self.has_cref:
            msg = 'Part has no reference curve.'
            raise ValueError(msg)
        return self.cref.p1

    @property
    def p2(self):
        """
        :return: The last point of the part reference curve.
        :rtype: afem.geometry.entities.Point

        :raises ValueError: If the part has no reference curve.
        """
        if not self.has_cref:
            msg = 'Part has no reference curve.'
            raise ValueError(msg)
        return self.cref.p2

    @property
    def nedges(self):
        """
        :return: The number of edges in the part shape.
        :rtype: int
        """
        return len(self.edges)

    @property
    def edges(self):
        """
        :return: All the edges of the part shape.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return ExploreShape.get_edges(self)

    @property
    def nfaces(self):
        """
        :return: The number of faces in the part shape.
        :rtype: int
        """
        return len(self.faces)

    @property
    def faces(self):
        """
        :return: All the faces of the part shape.
        :rtype: list[OCC.TopoDS.TopoDS_Face]
        """
        return ExploreShape.get_faces(self)

    def set_cref(self, cref):
        """
        Set the part reference curve.

        :param afem.geometry.entities.Curve cref: The curve.

        :return: None.

        :raise TypeError: If *cref* is an invalid curve.
        """
        if not CheckGeom.is_curve_like(cref):
            msg = 'Invalid curve type.'
            raise TypeError(msg)
        self._cref = cref

    def set_sref(self, sref):
        """
        Set the part reference surface.

        :param afem.geometry.entities.Surface sref: The surface.

        :return: None.

        :raise TypeError: If *sref* is an invalid surface.
        """
        if not CheckGeom.is_surface_like(sref):
            msg = 'Invalid surface type.'
            raise TypeError(msg)
        self._sref = sref

    def nullify(self):
        """
        Destroy reference of the part shape.

        :return: None
        """
        self.Nullify()

    def add_metadata(self, key, value):
        """
        Add metadata to the part.

        :param key: The key.
        :param value: The value.

        :return: None.
        """
        self._metadata[key] = value

    def get_metadata(self, key):
        """
        Get metadata.

        :param key: The key.

        :return: The key or *None* if not present.
        :rtype: object or None
        """
        try:
            return self._metadata[key]
        except KeyError:
            return None

    def add_subpart(self, key, subpart):
        """
        Add a sub-part to the part.

        :param str key: The key.
        :param afem.structure.entities.Part subpart: The sub-part.

        :return: None.

        :raise TypeError: If *subpart* is not a part.
        """
        if not isinstance(subpart, Part):
            msg = 'Sub-part is not a part.'
            raise TypeError(msg)
        self._subparts[key] = subpart

    def get_subpart(self, key):
        """
        Get a sub-part.

        :param str key: The key.

        :return: The sub-part. Returns *None* if the key is not present.
        :rtype: afem.structure.entities.Part or None
        """
        try:
            return self._subparts[key]
        except KeyError:
            return None

    def set_shape(self, shape):
        """
        Set the shape of the part.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.

        :return: None.

        :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
        """
        shape = CheckShape.to_shape(shape)
        if not shape:
            msg = 'Invalid shape.'
            raise RuntimeError(msg)

        self.TShape(shape.TShape())
        self.Location(shape.Location())
        self.Orientation(shape.Orientation())

    def local_to_global_u(self, u):
        """
        Convert local parameter from 0 <= u <= 1 to u1 <= u <= u2 using the
        part reference curve.

        :param float u: Local u-parameter.

        :return: Global u-parameter.
        :rtype: float

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)
        return self._cref.local_to_global_param(u)

    def eval_cref(self, u):
        """
        Evaluate point on reference curve.

        :param float u: The parameter.

        :return: The point.
        :rtype: afem.geometry.entities.Point

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)
        return self._cref.eval(u)

    def eval_sref(self, u, v):
        """
        Evaluate point on reference surface.

        :param float u: The u-parameter.
        :param float v: The v-parameter.

        :return: The point.
        :rtype: afem.geometry.entities.Point

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)
        return self._sref.eval(u, v)

    def point_from_parameter(self, ds, u0=None, is_local=False):
        """
        Evaluate point on reference curve at a distance from a parameter.

        :param float ds: The distance.
        :param float u0: The parameter. If not provided the first parameter
            of the reference curve will be used.
        :param bool is_local: Option specifying if the parameter is local or
            global.

        :return: The point.
        :rtype: afem.geometry.entities.Point

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)
        if u0 is None:
            u0 = self.cref.u1
        if is_local:
            u0 = self.local_to_global_u(u0)
        return PointFromParameter(self.cref, u0, ds).point

    def points_by_number(self, n, d1=None, d2=None, shape1=None,
                         shape2=None, tol=1.0e-7):
        """
        Create a specified number of points along the reference curve.

        :param int n: Number of points to create (*n* > 0).
        :param float d1: An offset distance for the first point. This is
            typically a positive number indicating a distance from *u1*
            towards *u2*.
        :param float d2: An offset distance for the last point. This is
            typically a negative number indicating a distance from *u2*
            towards *u1*.
        :param OCC.TopoDS.TopoDS_Shape shape1: A shape to define the first
            point. This shape is intersected with the edge or wire.
        :param OCC.TopoDS.TopoDS_Shape shape2: A shape to define the last
            point. This shape is intersected with the edge or wire.
        :param float tol: Tolerance.

        :return: The points.
        :rtype: list[afem.geometry.entities.Point]

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        edge = CheckShape.to_edge(self.cref)
        builder = PointsAlongShapeByNumber(edge, n, d1, d2, shape1, shape2,
                                           tol)
        return builder.points

    def points_by_distance(self, maxd, nmin=0, d1=None, d2=None, shape1=None,
                           shape2=None, tol=1.0e-7):
        """
        Create a points along the reference curve by distance.

        :param float maxd: The maximum allowed spacing between points. The
            actual spacing will be adjusted to not to exceed this value.
        :param int nmin: Minimum number of points to create.
        :param float d1: An offset distance for the first point. This is
            typically a positive number indicating a distance from *u1*
            towards *u2*.
        :param float d2: An offset distance for the last point. This is
            typically a negative number indicating a distance from *u2*
            towards *u1*.
        :param OCC.TopoDS.TopoDS_Shape shape1: A shape to define the first
            point. This shape is intersected with the edge or wire.
        :param OCC.TopoDS.TopoDS_Shape shape2: A shape to define the last
            point. This shape is intersected with the edge or wire.
        :param float tol: Tolerance.

        :return: The points.
        :rtype: list[afem.geometry.entities.Point]

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        edge = CheckShape.to_edge(self.cref)
        builder = PointsAlongShapeByDistance(edge, maxd, d1, d2, shape1,
                                             shape2, nmin, tol)
        return builder.points

    def point_to_cref(self, pnt, direction=None):
        """
        Project a point to reference curve.

        :param afem.geometry.entities.Point pnt: The point. Position will be
            updated.
        :param vector_like direction: Projection direction.

        :return: *True* if projected, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        proj = ProjectPointToCurve(pnt, self.cref, direction)
        if not proj.success:
            return False

        p = proj.nearest_point
        pnt.set_xyz(p.xyz)
        return True

    def points_to_cref(self, pnts, direction=None):
        """
        Project points to the reference curve.

        :param list[afem.geometry.entities.Point] pnts: The points. Position
            will be updated.
        :param vector_like direction: Projection direction.

        :return: List of status for each point.
        :rtype: list[bool]

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        success = []
        for p in pnts:
            status = self.point_to_cref(p, direction)
            success.append(status)
        return success

    def point_to_sref(self, pnt, direction=None):
        """
        Project a point to reference surface.

        :param afem.geometry.entities.Point pnt: The point. Position will be
            updated.
        :param vector_like direction: Projection direction.

        :return: *True* if projected, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        proj = ProjectPointToSurface(pnt, self.sref, direction)
        if not proj.success:
            return False

        p = proj.nearest_point
        pnt.set_xyz(p.xyz)
        return True

    def points_to_sref(self, pnts, direction=None):
        """
        Project points to reference surface.

        :param list[afem.geometry.entities.Point] pnts: The points. Position
            will be updated.
        :param vector_like direction: Projection direction.

        :return: List of status for each point.
        :rtype: list[bool]

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        success = []
        for p in pnts:
            status = self.point_to_sref(p, direction)
            success.append(status)
        return success

    def get_plane(self, u=None, ds=None, pnt=None, ref_pln=None):
        """
        Get a plane along the reference curve.

        :param u:
        :param ds:
        :param pnt:
        :param ref_pln:

        :return:

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        return CreateGeom.plane_on_curve(self.cref, u, ds, pnt, ref_pln)

    def invert_cref(self, pnt):
        """
        Invert the point on the reference curve.

        :param point_like pnt: The point.

        :return: The parameter on the reference curve. Returns *None* if no
            projection is found.
        :rtype: float or None

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        proj = ProjectPointToCurve(pnt, self.cref)
        if proj.success:
            return proj.nearest_param
        return None

    def invert_sref(self, pnt):
        """
        Invert the point on the reference surface.

        :param point_like pnt: The point.

        :return: The parameters on the reference surface (u, v). Returns
            *None* if no projection is found.
        :rtype: tuple(float) or tuple(None)

        :raise AttributeError: If part does not have a reference surface.
        """
        if not self.has_sref:
            msg = 'Part does not have a reference surface.'
            raise AttributeError(msg)

        proj = ProjectPointToSurface(pnt, self.sref)
        if proj.success:
            return proj.nearest_param
        return None, None

    def distance(self, other):
        """
        Find the minimum distance between the part and other shape.

        :param other: Other part or shape.
        :type other: OCC.TopoDS.TopoDS_Shape or afem.structure.entities.Part

        :return: The minimum distance and the points of minimum distance on
            each shape (dmin, p1, p2).
        :rtype: tuple(float, afem.geometry.entities.Point,
            afem.geometry.entities.Point)
        """
        # TODO Min distance for shapes.
        return ShapeTools.min_distance(self, other)

    def check(self):
        """
        Check the shape of the part.

        :return: *True* if shape is valid, *False* if not.
        :rtype: bool
        """
        return CheckShape.is_valid(self)

    def fix(self, min_tol=None, max_tol=None):
        """
        Attempt to fix the shape of the part.

        :param float min_tol: Minimum tolerance.
        :param float max_tol: Maximum tolerance.

        :return: None.
        """
        new_shape = FixShape(self, min_tol, max_tol).shape
        self.set_shape(new_shape)

    def reshape(self, tool):
        """
        Reshape the part shape with a tool.

        :param tool: The tool. It should provide the typical generated,
            modified, and deleted methods.

        :return: *True* if modified, *False* if not.
        :rtype: bool
        """
        return reshape_parts(tool, [self])

    def cut(self, cutter):
        """
        Cut the part shape.

        :param cutter: The cutter.
        :type cutter: OCC.TopoDS.TopoDS_Shape or afem.structure.entities.Part

        :return: *True* if shape was cut, *False* if not.
        :rtype: bool
        """
        cutter = CheckShape.to_shape(cutter)
        return cut_part(self, cutter)

    def split(self, splitter, split_both=True):
        """
        Split the part shape.

        :param splitter: The splitter.
        :type splitter: OCC.TopoDS.TopoDS_Shape or afem.structure.entities.Part
        :param bool split_both: Option to split both if *splitter* is a part.

        :return:
        """
        return split_part(self, splitter, split_both)


class CurvePart(Part):
    """
    Base class for curve parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None):
        super(CurvePart, self).__init__(label, shape, cref, None)

    @property
    def reshapes(self):
        """
        :return: The shapes to use when reshaping. For a curve part these
            are edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return self.edges

    @property
    def length(self):
        """
        :return: The length of all the edges of the part.
        :rtype: float
        """
        return LinearProps(self).length


class Beam(CurvePart):
    """
    Beam.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None):
        super(Beam, self).__init__(label, shape, cref)


class SurfacePart(Part):
    """
    Base class for surface parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(SurfacePart, self).__init__(label, shape, cref, sref)

    @property
    def reshapes(self):
        """
        :return: The shapes to use when reshaping. For a surface part these
            are faces.
        :rtype: list[OCC.TopoDS.TopoDS_Face]
        """
        return self.faces

    @property
    def area(self):
        """
        :return: The area of all faces of the part.
        :rtype: float
        """
        return SurfaceProps(self).area

    @property
    def stiffeners(self):
        """
        :return: List of stiffeners associated to the part.
        :rtype: list[afem.structure.entities.Stiffener1D]
        """
        return [part for part in self.subparts]

    @property
    def elements(self):
        """
        :return: The shell elements of the part.
        :rtype: set(afem.fem.entities.Elm2D)
        """
        smesh_mesh = MeshData.get_mesh().smesh_obj
        if not isinstance(smesh_mesh, SMESH_Mesh):
            return set()

        compound = ExploreShape.get_faces(self, True)
        submesh = smesh_mesh.GetSubMesh(compound)
        if submesh.IsEmpty():
            return set()

        submesh_ds = submesh.GetSubMeshDS()
        if not submesh_ds:
            return set()

        elm_iter = submesh_ds.GetElements()
        elm_set = set()
        while elm_iter.more():
            elm = Elm2D(elm_iter.next())
            elm_set.add(elm)

        return elm_set

    @property
    def nodes(self):
        """
        :return: The nodes of part.
        :rtype: set(afem.fem.entities.Node)
        """
        node_set = set()
        for e in self.elements:
            for n in e.nodes:
                node_set.add(n)
        return node_set

    def fuse(self, *other_parts):
        """
        Fuse with other parts.

        :param afem.structure.entities.Part other_parts: The other part(s)

        :return: *True* if fused, *False* if not.
        :rtype: bool
        """
        return fuse_surface_part(self, *other_parts)

    def sew(self, *other_parts):
        """
        Sew with other parts.

        :param afem.structure.entities.SurfacePart other_parts: The other
            part(s)

        :return: *True* if sewed, *False* if not.
        :rtype: bool
        """
        return sew_surface_parts([self] + list(other_parts))

    def merge(self, other, unify=False):
        """
        Merge other part or shape with this one.

        :param other: The other part or shape.
        :type other: afem.structure.entities.SurfacePart or
            OCC.TopoDS.TopoDS_Shape
        :param bool unify: Option to attempt to unify same domains.

        :return: *True* if merged, *False* if not.
        :rtype: bool
        """
        return merge_surface_part(self, other, unify)

    def unify(self, edges=True, faces=True, bsplines=False):
        """
        Attempt to unify the same domains of the part shape.

        :param bool edges: Option to unify all possible edges.
        :param bool faces: Option to unify all possible faces.
        :param bool bsplines: Option to concatenate the curves of edges if they
            are C1 continuous.

        :return: *True* if unified, *False* if not.
        :rtype: bool
        """
        return unify_surface_part(self, edges, faces, bsplines)

    def discard_by_solid(self, solid, tol=1.0e-7):
        """
        Discard faces of the part using a solid. Any faces of the part that
        have centroids inside the solid will be removed.

        :param OCC.TopoDS.TopoDS_Solid solid: The solid.
        :param float tol: The tolerance.

        :return: *True* if faces were discard, *False* if not.
        :rtype: bool
        """
        return discard_faces_by_solid(self, solid, tol)

    def discard_by_distance(self, shape, dmax):
        """
        Discard faces of the part using a shape and a distance. If the
        distance between a face of the part and the given shape is greater
        than the tolerance, then the face is removed.

        :param OCC.TopoDS.TopoDS_Shape shape: The solid.
        :param float dmax: The maximum distance allowed.

        :return: *True* if faces were discard, *False* if not.
        :rtype: bool
        """
        # TODO Support use of geometry.
        return discard_faces_by_distance(self, shape, dmax)

    def discard_by_cref(self):
        """
        Discard faces of the part by using the reference curve. An infinite
        solid is created at each end of the reference curve using the curve
        tangent. Any face that has a centroid in these solids is removed.

        :return: *True* if faces were discard, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        return discard_wing_part_faces(self)

    def shared_edges(self, other):
        """
        Get edges shared between the two parts.

        :param other: The other part or shape.
        :type other: afem.structure.entities.Part or OCC.TopoDS.TopoDS_Shape

        :return: Shared edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return get_shared_edges(self, other)

    def shared_nodes(self, other_part):
        """
        Get nodes shared between the two parts.

        :param afem.structure.entities.Part other_part: The other part.

        :return: Shared nodes.
        :rtype: list[afem.fem.entities.Node]
        """
        return get_shared_nodes(self, other_part)

    def cut_hole(self, ds, r):
        """
        Cut a circular hole in the part (in development).

        :param float ds: The distance along the reference curve from its
            first point.
        :param float r: The radius of the hole.

        :return: *True* if hole was cut, *False* if not.
        :rtype: bool

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        return cut_wing_part_with_circle(self, ds, r)


class WingPart(SurfacePart):
    """
    Base class for wing parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(WingPart, self).__init__(label, shape, cref, sref)


class Spar(WingPart):
    """
    Wing spar.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Spar, self).__init__(label, shape, cref, sref)


class Rib(WingPart):
    """
    Wing rib.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Rib, self).__init__(label, shape, cref, sref)


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(FuselagePart, self).__init__(label, shape, cref, sref)


class Bulkhead(FuselagePart):
    """
    Bulkhead.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Bulkhead, self).__init__(label, shape, cref, sref)


class Floor(FuselagePart):
    """
    Floor.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Floor, self).__init__(label, shape, cref, sref)


class Frame(FuselagePart):
    """
    Frame.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Frame, self).__init__(label, shape, cref, sref)


class Skin(SurfacePart):
    """
    Skin.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Skin, self).__init__(label, shape, cref, sref)


class Stiffener1D(CurvePart):
    """
    1-D stiffener for surface parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None):
        super(Stiffener1D, self).__init__(label, shape, cref)


class Stiffener2D(SurfacePart):
    """
    2-D stiffener for surface parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Stiffener2D, self).__init__(label, shape, cref, sref)


class Stringer(SurfacePart):
    """
    Stringer.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Stringer, self).__init__(label, shape, cref, sref)
