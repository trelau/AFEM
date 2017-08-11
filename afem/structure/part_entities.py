from OCC.SMESH import SMESH_Mesh
from OCC.TopoDS import TopoDS_Shape
from numpy import mean

from afem.config import logger
from afem.fem.elements import Elm2D
from afem.fem.meshes import MeshData
from afem.geometry.geom_check import CheckGeom
from afem.geometry.geom_create import (PlaneByNormal, PlaneFromParameter,
                                       PointFromParameter)
from afem.geometry.geom_project import (ProjectPointToCurve,
                                        ProjectPointToSurface)
from afem.graphics.graphics_viewer import ViewableItem
from afem.topology.shape_bop import CutShapes, FuseShapes, SplitShapes
from afem.topology.shape_check import CheckShape, ClassifyPointInSolid
from afem.topology.shape_create import (CompoundByShapes, FaceBySurface,
                                        HalfspaceByShape,
                                        PointsAlongShapeByDistance,
                                        PointsAlongShapeByNumber)
from afem.topology.shape_distance import DistanceShapeToShape
from afem.topology.shape_explore import ExploreShape
from afem.topology.shape_modify import (FixShape, RebuildShapeByTool,
                                        RebuildShapeWithShapes, SewShape,
                                        UnifyShape)
from afem.topology.shape_props import LinearProps, SurfaceProps

__all__ = ["Part", "CurvePart", "Beam", "SurfacePart", "WingPart", "Spar",
           "Rib", "FuselagePart", "Bulkhead", "Floor", "Frame", "Skin",
           "Stiffener1D", "Stiffener2D", "Stringer"]


class Part(TopoDS_Shape, ViewableItem):
    """
    Base class for all parts.

    :param str label: The label.
    :param OCC.TopoDS.TopoDS_Shape shape: The shape.
    :param cref: The reference curve.
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
        if cref is not None:
            self.set_cref(cref)
        if sref is not None:
            self.set_sref(sref)

        msg = ' '.join(['Creating part:', label])
        logger.info(msg)

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
        :rtype: list[afem.structure.part_entities.Part]
        """
        return self._subparts.values()

    @property
    def cref(self):
        """
        :return: The part reference curve.
        :rtype: afem.geometry.geom_entities.Curve
        """
        return self._cref

    @property
    def sref(self):
        """
        :return: The part reference surface.
        :rtype: afem.geometry.geom_entities.Surface
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
        :rtype: afem.geometry.geom_entities.Point

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
        :rtype: afem.geometry.geom_entities.Point

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

        :param afem.geometry.geom_entities.Curve cref: The curve.

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

        :param afem.geometry.geom_entities.Surface sref: The surface.

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
        :param afem.structure.part_entities.Part subpart: The sub-part.

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
        :rtype: afem.structure.part_entities.Part or None
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
        :rtype: afem.geometry.geom_entities.Point

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
        :rtype: afem.geometry.geom_entities.Point

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
        :rtype: afem.geometry.geom_entities.Point

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
        :rtype: list[afem.geometry.geom_entities.Point]

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
        :rtype: list[afem.geometry.geom_entities.Point]

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

        :param afem.geometry.geom_entities.Point pnt: The point. Position will be
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

        :param list[afem.geometry.geom_entities.Point] pnts: The points. Position
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

        :param afem.geometry.geom_entities.Point pnt: The point. Position will be
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

        :param list[afem.geometry.geom_entities.Point] pnts: The points. Position
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

    def get_plane(self, u0, ds, ref_pln=None, tol=1.0e-7):
        """
        Get a plane along the reference curve.

        :param float u0: The initial parameter.
        :param float ds: The distance along the curve from the given parameter.
        :param afem.geometry.geom_entities.Plane ref_pln: The normal of this plane
            will be used to define the normal of the new plane. If no plane is
            provided, then the first derivative of the curve will define the
            plane normal.
        :param float tol: Tolerance.

        :return: The plane.
        :rtype: afem.geometry.geom_entities.Plane

        :raise AttributeError: If part does not have a reference curve.
        """
        if not self.has_cref:
            msg = 'Part does not have a reference curve.'
            raise AttributeError(msg)

        return PlaneFromParameter(self.cref, u0, ds, ref_pln, tol)

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
        :type other: OCC.TopoDS.TopoDS_Shape or afem.structure.part_entities.Part

        :return: The minimum distance.
        :rtype: float
        """
        return DistanceShapeToShape(self, other).dmin

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

    def cut(self, cutter):
        """
        Cut the part shape and rebuild this part.

        :param cutter: The cutter.
        :type cutter: OCC.TopoDS.TopoDS_Shape or afem.structure.part_entities.Part

        :return: *True* if shape was cut, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        cut = CutShapes(self, cutter)
        if not cut.is_done:
            return False

        # Rebuild the parts
        if isinstance(self, CurvePart):
            rebuild = RebuildShapeByTool(self, cut, 'edge')
        elif isinstance(self, SurfacePart):
            rebuild = RebuildShapeByTool(self, cut, 'face')
        else:
            msg = 'Invalid part type in split operation.'
            raise TypeError(msg)

        new_shape = rebuild.new_shape
        self.set_shape(new_shape)

        return True

    def split(self, splitter, rebuild_both=True):
        """
        Split the part shape and rebuild this part. Optionally rebuild the
        splitter if it is a part. This method should handle splitting with
        parts and shapes or different types (i.e., splitting a face with an
        edge).

        :param splitter: The splitter.
        :type splitter: OCC.TopoDS.TopoDS_Shape or afem.structure.part_entities.Part
        :param bool rebuild_both: Option to rebuild both if *splitter* is a
            part.

        :return: *True* if split, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is or the splitter not a curve or
            surface part.
        """
        split = SplitShapes()
        split.add_arg(self)
        if rebuild_both:
            split.add_arg(splitter)
        else:
            split.add_tool(splitter)
        split.build()
        if not split.is_done:
            return False

        parts = [self]
        if rebuild_both:
            parts += [splitter]
        for part in parts:
            if isinstance(part, CurvePart):
                rebuild = RebuildShapeByTool(part, split, 'edge')
            elif isinstance(part, SurfacePart):
                rebuild = RebuildShapeByTool(part, split, 'face')
            else:
                msg = 'Invalid part type in split operation.'
                raise TypeError(msg)
            new_shape = rebuild.new_shape
            part.set_shape(new_shape)
        return True

    def rebuild(self, tool):
        """
        Rebuild the part shape with a supported tool.

        :param tool: The tool. It should provide the typical generated,
            modified, and deleted methods.

        :return: *True* if modified, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        if isinstance(self, CurvePart):
            rebuild = RebuildShapeByTool(self, tool, 'edge')
        elif isinstance(self, SurfaceProps):
            rebuild = RebuildShapeByTool(self, tool, 'face')
        else:
            msg = 'Invalid part type in split operation.'
            raise TypeError(msg)

        new_shape = rebuild.new_shape
        self.set_shape(new_shape)
        return True

    def discard_by_solid(self, solid, tol=None):
        """
        Discard shapes of the part using a solid. Any shapes of the part that
        have centroids inside the solid will be removed. Edges are checked
        for curve parts and faces are checked for surface parts.

        :param OCC.TopoDS.TopoDS_Solid solid: The solid.
        :param float tol: The tolerance. If not provided then the part
            tolerance will be used.

        :return: *True* if shapes were discarded, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        if isinstance(self, CurvePart):
            shapes = ExploreShape.get_edges(self)
        elif isinstance(self, SurfacePart):
            shapes = ExploreShape.get_faces(self)
        else:
            msg = 'Invalid part type in discard operation.'
            raise TypeError(msg)

        if tol is None:
            tol = self.tol

        rebuild = RebuildShapeWithShapes(self)
        classifer = ClassifyPointInSolid(solid, tol=tol)

        modified = False
        for shape in shapes:
            if isinstance(self, CurvePart):
                cg = LinearProps(shape).cg
            else:
                cg = SurfaceProps(shape).cg

            classifer.perform(cg, tol)
            if classifer.is_in:
                rebuild.remove(shape)
                modified = True

        if not modified:
            return False

        new_shape = rebuild.apply()
        self.set_shape(new_shape)
        return True

    def discard_by_distance(self, shape, dmax):
        """
        Discard shapes of the part using a shape and a distance. If the
        distance between a shape of the part and the given shape is greater
        than the tolerance, then the shape is removed. Edges are checked
        for curve parts and faces are checked for surface parts.

        :param OCC.TopoDS.TopoDS_Shape shape: The shape.
        :param float dmax: The maximum distance allowed.

        :return: *True* if shapes were discarded, *False* if not.
        :rtype: bool

        :raise TypeError: If this part is not a curve or surface part.
        """
        if isinstance(self, CurvePart):
            shapes = ExploreShape.get_edges(self)
        elif isinstance(self, SurfacePart):
            shapes = ExploreShape.get_faces(self)
        else:
            msg = 'Invalid part type in discard operation.'
            raise TypeError(msg)

        rebuild = RebuildShapeWithShapes(self)

        modified = False
        for part_shape in shapes:
            dmin = DistanceShapeToShape(shape, part_shape).dmin
            if dmin <= dmax:
                rebuild.remove(shape)
                modified = True

        if not modified:
            return False

        new_shape = rebuild.apply()
        self.set_shape(new_shape)
        return True


class CurvePart(Part):
    """
    Base class for curve parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param cref: The reference curve.
    :type cref: afem.geometry.geom_entities.Curve or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None):
        super(CurvePart, self).__init__(label, shape, cref, None)

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
    :type cref: afem.geometry.geom_entities.Curve or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

    :raise RuntimeError: If *shape* is not a valid shape or cannot be
            converted to a shape.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(SurfacePart, self).__init__(label, shape, cref, sref)

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
        :rtype: list[afem.structure.part_entities.Stiffener1D]
        """
        return [part for part in self.subparts]

    @property
    def elements(self):
        """
        :return: The shell elements of the part.
        :rtype: set(afem.fem.entities.Elm2D)
        """
        # TODO Support nodes and elements for both curve and surface parts.
        smesh_mesh = MeshData.get_mesh().smesh_obj
        if not isinstance(smesh_mesh, SMESH_Mesh):
            return set()

        faces = ExploreShape.get_faces(self)
        compound = CompoundByShapes(faces).compound
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
        Fuse with other surface parts and rebuild both.

        :param afem.structure.part_entities.SurfacePart other_parts: The other
            part(s).

        :return: *True* if fused, *False* if not.
        :rtype: bool
        """
        # Putting the other parts in a compound avoids fusing them to each
        # other.
        other_parts = list(other_parts)
        other_compound = CompoundByShapes(other_parts).compound

        fuse = FuseShapes(self, other_compound)
        if not fuse.is_done:
            return False

        # Rebuild the parts
        for part in [self] + other_parts:
            rebuild = RebuildShapeByTool(part, fuse)
            new_shape = rebuild.new_shape
            part.set_shape(new_shape)

        return True

    def sew(self, *other_parts):
        """
        Sew with other parts and rebuild all parts.

        :param afem.structure.part_entities.SurfacePart other_parts: The other
            part(s).

        :return: *True* if sewed, *False* if not.
        :rtype: bool
        """
        parts = [self] + list(other_parts)

        tol = mean([ExploreShape.get_tolerance(part, 0) for part in parts],
                   dtype=float)
        max_tol = max([ExploreShape.get_tolerance(part, 1) for part in parts])

        sew = SewShape(tol=tol, max_tol=max_tol, cut_free_edges=True,
                       non_manifold=True)
        for part in parts:
            sew.add(part)
        sew.perform()

        for part in parts:
            if not sew.is_modified(part):
                continue
            mod_shape = sew.modified(part)
            part.set_shape(mod_shape)
        return True

    def merge(self, other, unify=False):
        """
        Merge other surface part or shape with this one.

        :param other: The other part or shape.
        :type other: afem.structure.entities.SurfacePart or
            OCC.TopoDS.TopoDS_Shape
        :param bool unify: Option to attempt to unify same domains.

        :return: *True* if merged, *False* if not.
        :rtype: bool
        """
        # Fuse the parts
        fuse = FuseShapes(self, other)
        if not fuse.is_done:
            return False

        # Reset the shape
        self.set_shape(fuse.shape)

        # Unify if possible
        if not unify:
            return True
        return self.unify()

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
        unify = UnifyShape(self, edges, faces, bsplines)
        new_shape = unify.shape
        self.set_shape(new_shape)

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

        # Create vectors at each end of the reference curve pointing "out" of
        # the part
        u1, u2 = self.cref.u1, self.cref.u2
        v1 = self.cref.deriv(u1, 1)
        v2 = self.cref.deriv(u2, 1)
        # Reverse v1 so it's "out" of the part
        v1.reverse()

        # Create planes at each end
        p1 = self.cref.eval(u1)
        p2 = self.cref.eval(u2)
        pln1 = PlaneByNormal(p1, v1).plane
        pln2 = PlaneByNormal(p2, v2).plane

        # Translate points to define half space
        pref1 = p1 + 100. * v1.xyz
        pref2 = p2 + 100. * v2.xyz
        f1 = FaceBySurface(pln1).face
        f2 = FaceBySurface(pln2).face
        hs1 = HalfspaceByShape(f1, pref1).solid
        hs2 = HalfspaceByShape(f2, pref2).solid

        # Discard by solid
        status1 = self.discard_by_solid(hs1)
        status2 = self.discard_by_solid(hs2)

        return status1 or status2

    def shared_edges(self, other):
        """
        Get edges shared between the two parts.

        :param other: The other part or shape.
        :type other: afem.structure.part_entities.Part or OCC.TopoDS.TopoDS_Shape

        :return: Shared edges.
        :rtype: list[OCC.TopoDS.TopoDS_Edge]
        """
        return ExploreShape.get_shared_edges(self, other)

    def shared_nodes(self, other):
        """
        Get nodes shared between the two parts.

        :param afem.structure.part_entities.SurfacePart other: The other part.

        :return: Shared nodes.
        :rtype: list[afem.fem.entities.Node]
        """
        nodes1 = self.nodes
        nodes2 = other.nodes
        return list(nodes1 & nodes2)


class WingPart(SurfacePart):
    """
    Base class for wing parts.

    :param str label: The label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param cref: The reference curve.
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
    :type cref: afem.geometry.geom_entities.Curve or None

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
    :type cref: afem.geometry.geom_entities.Curve or None
    :param sref: The reference surface.
    :type sref: afem.geometry.geom_entities.Surface or None

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
