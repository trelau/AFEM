from OCC.Geom import Geom_Curve, Geom_Surface
from OCC.gp import gp_Pnt, gp_Vec
from numpy import ndarray

from afem.geometry.entities import *

__all__ = ["CheckGeom"]


class CheckGeom(object):
    """
    Geometry checker.
    """

    @staticmethod
    def is_geom(geom):
        """
        Check if the entity is geometry.

        :param geom: An entity.

        :return: *True* if entity is geometry, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geometry)

    @staticmethod
    def is_point_like(geom):
        """
        Check if the entity is point_like.

        :param geom: An entity.

        :return: *True* if the entity is point_like, *False* if not.
        :rtype: bool
        """
        if isinstance(geom, Point):
            return True
        if isinstance(geom, gp_Pnt):
            return True
        if isinstance(geom, (tuple, list, ndarray)):
            return len(geom) == 3
        return False

    @staticmethod
    def is_point(geom):
        """
        Check if the entity is a Point.

        :param geom: An entity.
        :return: *True* if the entity is a Point, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point)

    @staticmethod
    def is_point2d(geom):
        """
        Check if the entity is a Point2D.

        :param geom: An entity.
        :return: *True* if the entity is a Point2D, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point2D)

    @staticmethod
    def to_point(geom):
        """
        Convert entity to a Point if possible.

        :param geom: An entity.

        :return: The entity if already a Point, or a new Point if it is
            array_like. Returns *None* otherwise.
        :rtype: :class:`.Point` or None
        """
        if isinstance(geom, Point):
            return geom
        elif isinstance(geom, gp_Pnt):
            return Point(geom.XYZ())
        elif CheckGeom.is_point_like(geom):
            return Point(geom[0], geom[1], geom[2])
        return None

    @staticmethod
    def to_points(geoms):
        """
        Convert entities to Points if possible.

        :param list geoms: List of point_like entities.

        :return: List of Point instances.
        :rtype: list
        """
        return [CheckGeom.to_point(p) for p in geoms if
                CheckGeom.is_point_like(p)]

    @staticmethod
    def is_vector(geom):
        """
        Check if the entity is a Vector.

        :param geom: An entity.
        :return: *True* if the entity is a Vector, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Vector)

    @staticmethod
    def to_vector(geom):
        """
        Convert entity to a Vector if possible.

        :param geom: An entity.

        :return: The entity if already a Vector, or a new Vector if it is
            array_like. Returns *None* otherwise.
        :rtype: :class:`.Vector` or None
        """
        if isinstance(geom, Vector):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Vector(geom[0], geom[1], geom[2])
        elif isinstance(geom, Direction):
            return Vector(geom)
        elif isinstance(geom, gp_Vec):
            return Vector(geom.XYZ())
        return None

    @staticmethod
    def is_direction(geom):
        """
        Check if the entity is a Direction.

        :param geom: An entity.
        :return: *True* if the entity is a Direction, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Direction)

    @staticmethod
    def to_direction(geom):
        """
        Convert entity to a Direction if possible.

        :param geom: An entity.

        :return: The entity if already a Direction, or a new Direction if it is
            array_like. Returns *None* otherwise.
        :rtype: :class:`.Direction` or None
        """
        if isinstance(geom, Direction):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Direction(geom[0], geom[1], geom[2])
        elif isinstance(geom, Vector):
            return Direction(geom)
        return None

    @staticmethod
    def is_plane(geom):
        """
        Check if the entity is a Plane.

        :param geom: An entity.
        :return: *True* if the entity is a Plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Plane)

    @staticmethod
    def is_line(geom):
        """
        Check if the entity is a Line.

        :param geom: An entity.
        :return: *True* if the entity is a Line, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Line)

    @staticmethod
    def is_curve(geom):
        """
        Check if the entity is a curve.

        :param geom: An entity.
        :return: *True* if the entity is a curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, NurbsCurve)

    @staticmethod
    def is_curve_like(geom):
        """
        Check if the entity is a curve or line.

        :param geom: An entity.
        :return: *True* if the entity is a line or curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geom_Curve)

    @staticmethod
    def is_curve2d(geom):
        """
        Check if the entity is a 2-D curve.

        :param geom: An entity.
        :return: *True* if the entity is a 2-D curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, NurbsCurve2D)

    @staticmethod
    def is_surface(geom):
        """
        Check if the entity is a surface.

        :param geom: An entity.
        :return: *True* if the entity is a surface, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, NurbsSurface)

    @staticmethod
    def is_surface_like(geom):
        """
        Check if the entity is a surface or a plane.

        :param geom: An entity.
        :return: *True* if the entity is a surface or plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geom_Surface)

    @staticmethod
    def is_axis3(geom):
        """
        Check if the entity is an Axis3.

        :param geom: An entity.

        :return: *True* if Axis3, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Axis3)
