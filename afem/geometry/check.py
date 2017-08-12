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
        Check if the entity is a :class:`.Point`.

        :param geom: An entity.
        :return: *True* if the entity is a Point, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point)

    @staticmethod
    def is_point2d(geom):
        """
        Check if the entity is a :class:`.Point2D`.

        :param geom: An entity.
        :return: *True* if the entity is a Point2D, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point2D)

    @staticmethod
    def to_point(geom):
        """
        Convert entity to a :class:`.Point` if possible.

        :param geom: An entity.

        :return: The entity if already a Point, or a new Point if it is
            array_like. Returns *None* otherwise.
        :rtype: afem.geometry.entities.Point or None
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
        Convert entities to points if possible.

        :param list[point_like] geoms: List of entities.

        :return: List of points.
        :rtype: list[afem.geometry.entities.Point]
        """
        return [CheckGeom.to_point(p) for p in geoms if
                CheckGeom.is_point_like(p)]

    @staticmethod
    def is_vector(geom):
        """
        Check if the entity is a :class:`.Vector`.

        :param geom: An entity.
        :return: *True* if the entity is a Vector, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Vector)

    @staticmethod
    def to_vector(geom):
        """
        Convert entity to a Vector if possible.

        :param vector_like geom: An entity.

        :return: The entity if already a Vector, or a new Vector if it is
            array_like. Returns *None* otherwise.
        :rtype: afem.geometry.entities.Vector or None
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
        Check if the entity is a :class:`.Direction`.

        :param geom: An entity.
        :return: *True* if the entity is a Direction, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Direction)

    @staticmethod
    def to_direction(geom):
        """
        Convert entity to a Direction if possible.

        :param vector_like geom: An entity.

        :return: The entity if already a Direction, or a new Direction if it is
            array_like. Returns *None* otherwise.
        :rtype: afem.geometry.entities.Direction or None
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
        Check if the entity is a :class:`.Plane`.

        :param geom: An entity.
        :return: *True* if the entity is a Plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Plane)

    @staticmethod
    def is_line(geom):
        """
        Check if the entity is a :class:`.Line`.

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
        Check if the entity is an OpenCASCADE Geom_Curve.

        :param geom: An entity.
        :return: *True* if the entity is a Geom_Curve, *False* if not.
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
        Check if the entity is an OpenCASCADE Geom_Surface.

        :param geom: An entity.
        :return: *True* if the entity is a Geom_Surface, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geom_Surface)

    @staticmethod
    def is_axis3(geom):
        """
        Check if the entity is an :class:`.Axis3`.

        :param geom: An entity.

        :return: *True* if Axis3, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Axis3)

    @staticmethod
    def nearest_point(p, pnts):
        """
        Find the point nearest to a given point.

        :param point_like p: The point.
        :param  list[point_like] pnts: List of points.

        :return: The nearest point.
        :rtype: afem.geometry.entities.Point
        """
        p = CheckGeom.to_point(p)
        pnts = [CheckGeom.to_point(pi) for pi in pnts]

        dmin = p.distance(pnts[0])
        pmin = pnts[0]
        for pi in pnts[1:]:
            di = p.distance(pi)
            if di < dmin:
                dmin = di
                pmin = pi
        return pmin
