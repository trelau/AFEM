# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
from OCCT.gp import gp_Pnt, gp_Vec, gp_Pnt2d
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
    def to_point(geom):
        """
        Convert entity to a :class:`.Point` if possible.

        :param geom: An entity.

        :return: The entity if already a Point, or a new Point if it is
            point_like.
        :rtype: afem.geometry.entities.Point

        :raise TypeError: If entity cannot be converted to a Point.
        """
        if geom is None:
            return None

        if isinstance(geom, Point):
            return geom
        elif isinstance(geom, gp_Pnt):
            return Point(geom.XYZ())
        elif CheckGeom.is_point_like(geom):
            return Point(geom[0], geom[1], geom[2])
        else:
            raise TypeError('Cannot convert to Point.')

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
    def is_point2d_like(geom):
        """
        Check if the entity is point2d_like.

        :param geom: An entity.

        :return: *True* if the entity is point2d_like, *False* if not.
        :rtype: bool
        """
        if isinstance(geom, Point2D):
            return True
        if isinstance(geom, gp_Pnt2d):
            return True
        if isinstance(geom, (tuple, list, ndarray)):
            return len(geom) == 2
        return False

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
    def to_point2d(geom):
        """
        Convert entity to a :class:`.Point2D` if possible.

        :param geom: An entity.

        :return: The entity if already a Point2D, or a new Point2D if it is
            point2d_like.
        :rtype: afem.geometry.entities.Point2D

        :raise TypeError: If entity cannot be converted to a Point2D.
        """
        if geom is None:
            return None

        if isinstance(geom, Point2D):
            return geom
        elif isinstance(geom, gp_Pnt2d):
            return Point2D(geom.XY())
        elif CheckGeom.is_point2d_like(geom):
            return Point2D(geom[0], geom[1])
        else:
            raise TypeError('Cannot convert to Point2D.')

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
            vector_like.
        :rtype: afem.geometry.entities.Vector

        :raise TypeError: If entity cannot be converted to a Vector.
        """
        if geom is None:
            return None

        if isinstance(geom, Vector):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Vector(geom[0], geom[1], geom[2])
        elif isinstance(geom, Direction):
            return Vector(geom)
        elif isinstance(geom, gp_Vec):
            return Vector(geom.XYZ())
        else:
            raise TypeError('Cannot convert to Vector.')

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
            vector_like.
        :rtype: afem.geometry.entities.Direction

        :raise TypeError: If entity cannot be converted to a Direction.
        """
        if geom is None:
            return None

        if isinstance(geom, Direction):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Direction(geom[0], geom[1], geom[2])
        elif isinstance(geom, Vector):
            return Direction(geom)
        else:
            raise TypeError('Cannot convert to Direction.')

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
        return isinstance(geom, Curve)

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
        return isinstance(geom, Surface)

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
