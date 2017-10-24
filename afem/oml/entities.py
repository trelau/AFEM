from OCC.TopoDS import TopoDS_Solid

from afem.geometry.create import PlaneByPoints, TrimmedCurveByParameters
from afem.geometry.entities import NurbsCurve, Surface
from afem.geometry.project import (ProjectPointToCurve,
                                   ProjectPointToSurface)
from afem.graphics.viewer import ViewableItem
from afem.topology import ExploreShape
from afem.topology.bop import IntersectShapes
from afem.topology.check import CheckShape
from afem.topology.create import FaceBySurface, WiresByConnectedEdges
from afem.topology.distance import DistancePointToShapes
from afem.topology.entities import BBox
from afem.topology.modify import DivideC0Shape, DivideClosedShape

__all__ = ["Body"]


class Body(TopoDS_Solid, ViewableItem):
    """
    Base class for OML bodies.

    :param solid: The solid. It will be downcasted to a solid if a generic
        shape is provided.
    :type solid: OCC.TopoDS.TopoDS_Shape or OCC.TopoDS.TopoDS_Solid
    :param str name: The name.

    :raise TypeError: If *shape* is not a solid.
    """

    def __init__(self, solid, name=None):
        super(Body, self).__init__()
        ViewableItem.__init__(self)
        self.set_solid(solid)
        self._metadata = {}
        self._name = name
        self._sref = None
        self._sref_shape = None

    @property
    def name(self):
        """
        :return: The name.
        :rtype: str
        """
        return self._name

    @property
    def solid(self):
        """
        :return: The shape.
        :rtype: OCC.TopoDS.TopoDS_Solid
        """
        return self

    @property
    def outer_shell(self):
        """
        :return: The outer shell.
        :rtype: OCC.TopoDS.TopoDS_Shell
        """
        return ExploreShape.outer_shell(self)

    @property
    def metadata(self):
        """
        :return: The metadata dictionary.
        :rtype: dict
        """
        return self._metadata

    @property
    def sref(self):
        """
        :return: The reference surface.
        :rtype: afem.geometry.entities.Surface
        """
        return self._sref

    @property
    def sref_shape(self):
        """
        :return: The reference shape.
        :rtype: OCC.TopoDS.TopoDS_Shape
        """
        return self._sref_shape

    @property
    def u1(self):
        """
        :return: Lower u-parameter of reference surface.
        :rtype: float
        """
        return self.sref.u1

    @property
    def u2(self):
        """
        :return: Upper u-parameter of reference surface.
        :rtype: float
        """
        return self.sref.u2

    @property
    def v1(self):
        """
        :return: Lower v-parameter of reference surface.
        :rtype: float
        """
        return self.sref.v1

    @property
    def v2(self):
        """
        :return: Upper v-parameter of reference surface.
        :rtype: float
        """
        return self.sref.v2

    def set_name(self, name):
        """
        Set name of body.

        :param str name: The name.

        :return: None.
        """
        self._name = name

    def add_metadata(self, key, value):
        """
        Add metadata to the body.

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

    def set_solid(self, solid):
        """
        Set the shape for the OML body.

        :param OCC.TopoDS.TopoDS_Solid solid: The solid.

        :return: None.

        :raise TypeError: If *shape* is not a solid.
        """
        if not CheckShape.is_solid(solid):
            msg = 'Invalid shape provided to Body. Requires a TopoDS_Solid.'
            raise RuntimeError(msg)
        solid = CheckShape.to_solid(solid)
        self.TShape(solid.TShape())
        self.Location(solid.Location())
        self.Orientation(solid.Orientation())
        return True

    def set_sref(self, srf, divide_closed=True, divide_c0=True):
        """
        Set the reference surface.

        :param afem.geometry.entities.Surface srf: The surface.
        :param bool divide_closed: Option to divide the surface if closed
            when creating the reference shape.
        :param bool divide_c0: Option to divide the surface at C0 boundaries
            when creating the reference shape.

        :return: None.

        :raise TypeError: If *srf* is not a supported type.
        """
        if not isinstance(srf, Surface):
            msg = 'Unsupported surface type.'
            raise TypeError(msg)

        self._sref = srf

        shape = FaceBySurface(srf).face
        if divide_closed:
            shape = DivideClosedShape(shape).shape
        if divide_c0:
            shape = DivideC0Shape(shape).shape
        self._sref_shape = shape

    def eval(self, u, v):
        """
        Evaluate a point on the reference surface.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.

        :return: Point on reference surface.
        :rtype: afem.geometry.entities.Point
        """
        return self.sref.eval(u, v)

    def norm(self, u, v):
        """
        Evaluate the surface normal of the reference surface.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.

        :return: Reference surface normal.
        :rtype: afem.geometry.entities.Vector
        """
        return self.sref.norm(u, v)

    def invert(self, p):
        """
        Find the parameters on the reference surface by inverting the point.

        :param point_like p: The point.

        :return: Parameters on the wing reference surface (u, v).
        :rtype: tuple(float)

        :raise RuntimeError: If no points are found in the projection
            algorithm.
        """
        proj = ProjectPointToSurface(p, self.sref)
        if not proj.success:
            msg = 'Failed to invert point.'
            raise RuntimeError(msg)
        return proj.nearest_param

    def extract_plane(self, u1, v1, u2, v2):
        """
        Extract a plane between parameters on the reference surface. The
        plane will be defined by three points. The first point is at (u1, v1),
        the second point is at (u2, v2), and the third point will be offset
        from the first point in the direction of the reference surface
        normal. The points should not be collinear.

        :param float u1: First u-parameter.
        :param float v1: First v-parameter.
        :param float u2: Second u-parameter.
        :param float v2: Second v-parameter.

        :return: The plane.
        :rtype: afem.geometry.entities.Plane
        """
        p1 = self.eval(u1, v1)
        p2 = self.eval(u2, v2)
        vn = self.norm(u1, v1)
        p3 = p1.copy()
        p3.translate(vn)
        return PlaneByPoints(p1, p2, p3).plane

    def extract_curve(self, u1, v1, u2, v2, basis_shape=None,
                      refine_edges=True):
        """
        Extract a trimmed curve within the reference surface between the
        parameters.

        :param float u1: First u-parameter.
        :param float v1: First v-parameter.
        :param float u2: Second u-parameter.
        :param float v2: Second v-parameter.
        :param basis_shape: The shape that will be used to intersect with
            the reference shape. If not provided a plane will be
            created using the *extract_plane()* method. The parameters
            should create points that are on or very near the intersection
            between these two shapes. If they are not they will be projected to
            the intersection which could yield unanticipated results.
        :type basis_shape: afem.geometry.entities.Surface or
            OCC.TopoDS.TopoDS_Shape
        :param bool refine_edges: Option to refine the edges after the
            Boolean operation between the reference shape and the basis
            shape.

        :return: The curve.
        :rtype: afem.geometry.entities.TrimmedCurve

        :raise RuntimeError: If method fails.
        """
        p1 = self.eval(u1, v1)
        p2 = self.eval(u2, v2)

        if basis_shape is None:
            basis_shape = self.extract_plane(u1, v1, u2, v2)
        basis_shape = CheckShape.to_shape(basis_shape)

        bop = IntersectShapes(basis_shape, self.sref_shape)
        if refine_edges:
            bop.refine_edges()
        shape = bop.shape

        edges = ExploreShape.get_edges(shape, unique=True)
        builder = WiresByConnectedEdges(edges)
        if builder.nwires == 0:
            msg = 'Failed to extract any curves.'
            raise RuntimeError(msg)

        if builder.nwires == 1:
            wire = builder.wires[0]
        else:
            dist = DistancePointToShapes(p1, builder.wires)
            wire = dist.nearest_shape
        crv = ExploreShape.curve_of_shape(wire)
        if not isinstance(crv, NurbsCurve):
            msg = 'Unsupported curve type created.'
            raise RuntimeError(msg)

        proj = ProjectPointToCurve(p1, crv)
        if not proj.success:
            msg = 'Failed to project point to reference curve.'
            raise RuntimeError(msg)
        u1c = proj.nearest_param

        proj = ProjectPointToCurve(p2, crv)
        if not proj.success:
            msg = 'Failed to project point to reference curve.'
            raise RuntimeError(msg)
        u2c = proj.nearest_param

        if u1c > u2c:
            crv.reverse()
            u1c, u2c = crv.reversed_u(u1c), crv.reversed_u(u2c)

        return TrimmedCurveByParameters(crv, u1c, u2c).curve

    def isocurve(self, u=None, v=None):
        """
        Extract iso-curve in the reference surface.

        :param float u: Constant u-parameter. Will override *v* if both are
            provided.
        :param float v: Constant v-parameter.

        :return: The iso-curve.
        :rtype: afem.geometry.entities.Curve

        :raise TypeError: If both *u* and *v* are None.
        """
        if u is not None:
            return self.sref.u_iso(u)
        elif v is not None:
            return self.sref.v_iso(v)
        else:
            msg = 'Invalid parameter input.'
            raise TypeError(msg)

    def bbox(self, tol=None):
        """
        Return a bounding box of the body.

        :param tol: Optional tolerance to enlarge the bounding box.
        :type tol: float or None

        :return: Bounding box of the body.
        :rtype: afem.topology.entities.BBox
        """
        bbox = BBox()
        bbox.add_shape(self)
        if tol is not None:
            bbox.enlarge(tol)
        return bbox
