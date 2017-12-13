#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

from warnings import warn

from OCCT.BRepExtrema import BRepExtrema_DistShapeShape

from afem.geometry.check import CheckGeom
from afem.topology.check import CheckShape

__all__ = ["DistanceShapeToShape", "DistanceShapeToShapes",
           "DistancePointToShapes"]


class DistanceShapeToShape(object):
    """
    Calculate minimum distance between two shapes.

    :param OCCT.TopoDS.TopoDS_Shape shape1: The first shape.
    :param OCCT.TopoDS.TopoDS_Shape shape2: The second shape.

    :raise RuntimeError: If OCC method fails.

    Usage:

    >>> from afem.topology import *
    >>> v1 = VertexByPoint((0., 0., 0.)).vertex
    >>> v2 = VertexByPoint((10., 0., 0.)).vertex
    >>> tool = DistanceShapeToShape(v1, v2)
    >>> tool.nsol
    1
    >>> tool.dmin
    10.0
    """

    def __init__(self, shape1, shape2):
        self._tool = BRepExtrema_DistShapeShape(shape1, shape2)
        if not self._tool.IsDone():
            msg = 'OCC BRepExtrema_DistShapeShape failed.'
            raise RuntimeError(msg)

    @property
    def nsol(self):
        """
        :return: The number of solutions satisfying the minimum distance.
        :rtype:
        """
        return self._tool.NbSolution()

    @property
    def dmin(self):
        """
        :return: The minimum distance.
        :rtype: float
        """
        return self._tool.Value()


class DistanceShapeToShapes(object):
    """
    Calculate the minimum distance between a shape and other shapes. Sort the
    results by distance.

    :param OCCT.TopoDS.TopoDS_Shape shape: The main shape.
    :param list[OCCT.TopoDS.TopoDS_Shape] other_shapes: The other shapes.

    :raises RuntimeWarning: If the distance between two shapes cannot be
        found. This shape will be ignored and process will continue.

    Usage:

    >>> from afem.topology import *
    >>> v1 = VertexByPoint((0., 0., 0.)).vertex
    >>> v2 = VertexByPoint((5., 0., 0.)).vertex
    >>> v3 = VertexByPoint((10., 0., 0.)).vertex
    >>> tool = DistanceShapeToShapes(v1, [v3, v2])
    >>> tool.dmin
    5.0
    >>> tool.dmax
    10.0
    >>> tool.sorted_distances
    [5.0, 10.0]
    """

    def __init__(self, shape, other_shapes):
        results = []
        for shape2 in other_shapes:
            dist = DistanceShapeToShape(shape, shape2)
            if dist.nsol == 0:
                warn("Could not calculate distance to a shape. Continuing...",
                     RuntimeWarning)
                continue
            results.append((dist.dmin, shape2))

        results.sort(key=lambda tup: tup[0])
        self._distances = [data[0] for data in results]
        self._shapes = [data[1] for data in results]

    @property
    def dmin(self):
        """
        :return: The minimum distance of all shapes.
        :rtype: float
        """
        return self._distances[0]

    @property
    def dmax(self):
        """
        :return: The maximum distance of all shapes.
        :rtype: float
        """
        return self._distances[-1]

    @property
    def sorted_distances(self):
        """
        :return: List of sorted distances.
        :rtype: list[float]
        """
        return self._distances

    @property
    def nearest_shape(self):
        """
        :return: The nearest shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shapes[0]

    @property
    def farthest_shape(self):
        """
        :return: The farthest shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shapes[-1]

    @property
    def sorted_shapes(self):
        """
        :return: List of shapes sorted by distance.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return self._shapes


class DistancePointToShapes(DistanceShapeToShapes):
    """
    Calculate the minimum distance between a point and other shapes. Sort the
    results by distance. This method converts the point to a vertex and then
    uses :class:`.DistanceShapeToShapes`.

    :param point_like pnt: The point.
    :param list[OCCT.TopoDS.TopoDS_Shape] other_shapes: The other shapes.

    :raise TypeError: If *pnt* cannot be converted to a point.
    """

    def __init__(self, pnt, other_shapes):
        pnt = CheckGeom.to_point(pnt)
        if not pnt:
            msg = 'Invalid point type provided.'
            raise TypeError(msg)

        v = CheckShape.to_vertex(pnt)
        super(DistancePointToShapes, self).__init__(v, other_shapes)
