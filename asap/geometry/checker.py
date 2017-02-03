from OCC.Geom import Geom_Curve, Geom_Surface
from OCC.gp import gp_Pnt
from numpy import ndarray

from .axes import Axis3
from .curves import Line, NurbsCurve, NurbsCurve2D
from .geom import Geometry
from .points import Point, Point2D
from .surfaces import NurbsSurface, Plane
from .vectors import Direction, Vector


class CheckGeom(object):
    """
    Geometry checker.
    """

    @staticmethod
    def is_geom(geom):
        """
        Check if the entity is geometry.

        :param geom: Entity.

        :return: *True* if entity is geometry, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geometry)

    @staticmethod
    def is_point_like(geom):
        """
        Check to see if the entity is point_like, which includes a Point
        instance, a tuple, list, or NumPy array of shape 1 x 3 that contains
        pnt coordinates.

        :param geom: Entity
        :type geom: :class:`.Point`  or array_like

        :return: *True* if the entity is point_like, *False* if not.
        :rtype: bool
        """
        if isinstance(geom, Point):
            return True
        if isinstance(geom, (tuple, list, ndarray)):
            return len(geom) == 3
        return False

    @staticmethod
    def is_point(geom):
        """
        Check to see if the entity is a point.

        :param geom: Geometric entity.
        :return: *True* if the entity is a point, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point)

    @staticmethod
    def is_point2d(geom):
        """
        Check to see if the entity is a Point2D.

        :param geom: Geometric entity.
        :return: *True* if the entity is a Point2D, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point2D)

    @staticmethod
    def to_point(geom):
        """
        Check to see if the entity is a :class:`.Point` instance and return
        the instance if it is. If the entity is point_like, create a new
        :class:`.Point` instance and return it.

        :param geom: Geometric entity or possible array.

        :return: The Point instance if already a point, or a new Point in
            case it is array_like.
        :rtype: :class:`.Point`
        """
        if isinstance(geom, Point):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Point(*geom)
        elif isinstance(geom, gp_Pnt):
            return Point(geom.XYZ())
        return None

    @staticmethod
    def to_points(geoms):
        """
        Check to see if the entities are a :class:`.Point` instance or
        convert them if they are point_like.

        :param list geoms: List of point_like entities.

        :return: List of :class:`.Point` instances.
        :rtype: list
        """
        return [CheckGeom.to_point(p) for p in geoms]

    @staticmethod
    def is_vector(geom):
        """
        Check to see if the entity is a vector.

        :param geom: Geometric entity.
        :return: *True* if the entity is a vector, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Vector)

    @staticmethod
    def to_vector(geom):
        """
        Check to see if the entity is a :class:`.Vector` instance and return
        the instance if it is. If the entity is array_like, create a new
        :class:`.Vector` instance and return it.

        :param geom: Geometric entity or possible array.

        :return: The Vector.
        :rtype: :class:`.Vector`
        """
        if isinstance(geom, Vector):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Vector(*geom)
        elif isinstance(geom, Direction):
            return Vector(Direction)
        return None

    @staticmethod
    def is_direction(geom):
        """
        Check to see if the entity is a direction.

        :param geom: Geometric entity.
        :return: *True* if the entity is a direction, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Direction)

    @staticmethod
    def to_direction(geom):
        """
        Check to see if the entity is a :class:`.Direction` instance and return
        the instance if it is. If the entity is array_like, create a new
        :class:`.Direction` instance and return it.

        :param geom: Geometric entity or possible array.

        :return: The Direction.
        :rtype: :class:`.Direction`
        """
        if isinstance(geom, Direction):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Direction(*geom)
        elif isinstance(geom, Vector):
            return Direction(Vector)
        return None

    @staticmethod
    def is_plane(geom):
        """
        Check to see if the entity is a plane.

        :param geom: Geometric entity.
        :return: *True* if the entity is a plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Plane)

    @staticmethod
    def is_line(geom):
        """
        Check to see if the entity is a line.

        :param geom: Geometric entity.
        :return: *True* if the entity is a line, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Line)

    @staticmethod
    def is_curve(geom):
        """
        Check to see if the entity is a curve.

        :param geom: Geometric entity.
        :return: *True* if the entity is a curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, NurbsCurve)

    @staticmethod
    def is_curve_like(geom):
        """
        Check to see if the entity is a curve or line.

        :param geom: Geometric entity.
        :return: *True* if the entity is a line or curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geom_Curve)

    @staticmethod
    def is_curve2d(geom):
        """
        Check to see if the entity is a 2-D curve.

        :param geom: Geometric entity.
        :return: *True* if the entity is a curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, NurbsCurve2D)

    @staticmethod
    def is_surface(geom):
        """
        Check to see if the entity is a surface.

        :param geom: Geometric entity.
        :return: *True* if the entity is a surface, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, NurbsSurface)

    @staticmethod
    def is_surface_like(geom):
        """
        Check to see if the entity is a surface or a plane.

        :param geom: Geometric entity.
        :return: *True* if the entity is a surface or plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geom_Surface)

    @staticmethod
    def is_axis3(geom):
        """
        Check to see if the entity is an axis.

        :param geom:

        :return: *True* if axis, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Axis3)
