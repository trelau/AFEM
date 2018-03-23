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
from math import radians, tan
from warnings import warn

from OCCT.BRepAdaptor import BRepAdaptor_CompCurve

from afem.geometry.check import CheckGeom
from afem.geometry.create import *
from afem.geometry.entities import *
from afem.structure.entities import *
from afem.topology.bop import *
from afem.topology.check import CheckShape
from afem.topology.create import *
from afem.topology.distance import DistanceShapeToShape
from afem.topology.explore import ExploreFreeEdges, ExploreShape
from afem.topology.modify import SewShape
from afem.topology.offset import SweepShapeWithNormal, SweepShape
from afem.topology.props import *

__all__ = ["CurvePartByShape", "Beam1DByShape", "Beam1DByCurve",
           "Beam1DByPoints",
           "SurfacePartByShape", "SparByParameters", "SparByPoints",
           "SparByEnds",
           "SparBySurface", "SparByShape",
           "SparBetweenShapes", "SparsBetweenPlanesByNumber",
           "SparsBetweenPlanesByDistance", "SparsAlongCurveByNumber",
           "SparsAlongCurveByDistance",
           "RibByParameters",
           "RibByPoints", "RibBySurface", "RibByShape", "RibBetweenShapes",
           "RibByOrientation",
           "RibsBetweenPlanesByNumber", "RibsBetweenPlanesByDistance",
           "RibsAlongCurveByNumber", "RibsAlongCurveByDistance",
           "RibsAlongCurveAndSurfaceByDistance",
           "BulkheadBySurface", "FloorBySurface",
           "FrameByPlane", "FramesByPlanes", "FramesBetweenPlanesByNumber",
           "FramesBetweenPlanesByDistance", "SkinBySolid",
           "SkinByBody", "StringerByShape", "Beam2DBySweep",
           "CreatePartByName"]

_type_to_part = {
    'Part': Part,
    'CurvePart': CurvePart,
    'Beam1D': Beam1D,
    'SurfacePart': SurfacePart,
    'WingPart': WingPart,
    'Spar': Spar,
    'Rib': Rib,
    'FuselagePart': FuselagePart,
    'Bulkhead': Bulkhead,
    'Floor': Floor,
    'Frame': Frame,
    'Skin': Skin,
    'Stiffener1D': Stiffener1D,
    'Stiffener2D': Stiffener2D,
    'Stringer': Stringer,
    'Beam2D': Beam2D
}


class CreatePartByName(object):
    """
    Create a part by its type name. This is a tool mostly for loading models
    and it not recommended for general use.

    :param str type_: The type of part.
    """

    def __init__(self, type_, *args, **kwargs):
        self._part = _type_to_part[type_](*args, **kwargs)

    @property
    def part(self):
        """
        :return: The part.
        :rtype: afem.structure.entities.Part
        """
        return self._part


# CURVE PART ------------------------------------------------------------------

class CurvePartByShape(object):
    """
    Create a curve part using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: afem.geometry.entities.Curve or OCCT.TopoDS.TopoDS_Shape
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> from afem.structure import CurvePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> part = CurvePartByShape('part', e).curve_part
    """

    def __init__(self, label, shape, cref=None, group=None):
        if isinstance(shape, Curve):
            cref = shape
            shape = EdgeByCurve(shape).edge

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._curve_part = CurvePart(label, shape, cref, None, group)

    @property
    def curve_part(self):
        """
        :return: The curve part.
        :rtype: afem.structure.entities.CurvePart
        """
        return self._curve_part


# BEAM1D ----------------------------------------------------------------------

class Beam1DByShape(object):
    """
    Create a beam using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: afem.geometry.entities.Curve or OCCT.TopoDS.TopoDS_Shape
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> from afem.structure import Beam1DByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> beam = Beam1DByShape('part', e).beam
    """

    def __init__(self, label, shape, cref=None, group=None):
        if isinstance(shape, Curve):
            shape = EdgeByCurve(shape).edge

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._beam = Beam1D(label, shape, cref, None, group)

    @property
    def beam(self):
        """
        :return: The beam.
        :rtype: afem.structure.entities.Beam1D
        """
        return self._beam


class Beam1DByCurve(Beam1DByShape):
    """
    Create a beam using a curve.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.structure import Beam1DByCurve
    >>> c = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> part = Beam1DByCurve('part', c).beam
    """

    def __init__(self, label, crv, group=None):
        e = EdgeByCurve(crv).edge
        super(Beam1DByCurve, self).__init__(label, e, crv, group)


class Beam1DByPoints(Beam1DByCurve):
    """
    Create a beam between two points.

    :param str label: Part label.
    :param point_like p1: First point.
    :param point_like p2: Second point.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    Usage:

    >>> from afem.structure import Beam1DByPoints
    >>> beam = Beam1DByPoints('part', (0., 0., 0.), (10., 0., 0.)).beam
    >>> beam.length
    10.0
    """

    def __init__(self, label, p1, p2, group=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)
        c = NurbsCurveByPoints([p1, p2]).curve
        super(Beam1DByPoints, self).__init__(label, c, group)


# SURFACE PART ----------------------------------------------------------------

class SurfacePartByShape(object):
    """
    Create a surface part using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: afem.geometry.entities.Surface or OCCT.TopoDS.TopoDS_Shape
    :param afem.geometry.entities.Curve cref: The reference curve.
    :param afem.geometry.entities.Surface sref: The reference surface. If
        not provided then the basis surface of the largest face of the shape
        will be used.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    Usage:

    >>> from afem.topology import EdgeByPoints, FaceByDrag
    >>> from afem.structure import SurfacePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> f = FaceByDrag(e, (0., 10., 0.)).face
    >>> part = SurfacePartByShape('part', f).surface_part
    >>> part.area
    99.99999999999997
    """

    def __init__(self, label, shape, cref=None, sref=None, group=None):
        if isinstance(shape, Surface):
            sref = shape
            shape = FaceBySurface(shape).face

        if sref is None:
            sref = ExploreShape.surface_of_shape(shape)

        self._surface_part = SurfacePart(label, shape, cref, sref, group)

    @property
    def surface_part(self):
        """
        :return: The surface part.
        :rtype: afem.structure.entities.SurfacePart
        """
        return self._surface_part


# SPAR ------------------------------------------------------------------------

class SparByParameters(object):
    """
    Create a spar between wing parameters.

    :param str label: Part label.
    :param float u1: Starting point u-parameter.
    :param float v1: Starting point v-parameter.
    :param float u2: Ending point u-parameter.
    :param float v2: Ending point v-parameter.
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.

    .. note::

        * If a basis shape is provided, then the starting and ending parameters
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, u1, v1, u2, v2, body, basis_shape=None,
                 group=None):
        # Determine reference surface and basis shape
        if basis_shape is None:
            pln = body.extract_plane(u1, v1, u2, v2)
            sref = pln
            basis_shape = FaceBySurface(pln).face
        elif isinstance(basis_shape, Surface):
            sref = basis_shape
            basis_shape = FaceBySurface(sref).face
        else:
            sref = ExploreShape.surface_of_shape(basis_shape)

        # Extract cref
        cref = body.extract_curve(u1, v1, u2, v2, basis_shape)

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._spar = Spar(label, shape, cref, sref, group)

    @property
    def spar(self):
        """
        :return: The spar.
        :rtype: afem.structure.entities.Spar
        """
        return self._spar


class SparByPoints(SparByParameters):
    """
    Create a spar between two points. This method inverts the starting and
    ending points and then uses :class:`.SparByParameters`.

    :param str label: Part label.
    :param point_like p1: Starting point.
    :param point_like p2: End point.
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        * The starting and ending points should be on or near the body
          reference surface since they are projected to find parameters.

        * If a basis shape is provided, then the starting and ending points
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, p1, p2, body, basis_shape=None, group=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        # Invert points
        u1, v1 = body.invert(p1)
        u2, v2 = body.invert(p2)

        # Use SparByParameters
        super(SparByPoints, self).__init__(label, u1, v1, u2, v2, body,
                                           basis_shape, group)


class SparByEnds(SparByParameters):
    """
    Create a spar by defining its endpoints which can be either points or
    parameters on a reference surface.

    :param str label: Part label.
    :param e1: Start point as a point or surface parameters (u1, v1).
    :type e1: point_like or collections.Sequence[float]
    :param e2: End point as a point or surface parameters (u2, v2).
    :type e2: point_like or collections.Sequence[float]
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    :raise TypeError: If *e1* or *e2* are not *point_like* or a sequence of
        two surface parameters.
    """

    def __init__(self, label, e1, e2, body, basis_shape=None, group=None):
        if len(e1) == 2:
            u1, v1 = e1
        elif CheckGeom.is_point_like(e1):
            u1, v1 = body.invert(e1)
        else:
            raise TypeError('Invalid type for e1.')

        if len(e2) == 2:
            u2, v2 = e2
        elif CheckGeom.is_point_like(e2):
            u2, v2 = body.invert(e2)
        else:
            raise TypeError('Invalid type for e2.')

        super(SparByEnds, self).__init__(label, u1, v1, u2, v2, body,
                                         basis_shape, group)


class SparBySurface(object):
    """
    Create a spar using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *srf* and the body
        reference shape.
    """

    def __init__(self, label, srf, body, group=None):
        basis_shape = FaceBySurface(srf).face

        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape,
                                  approximate=True)
        edges = ExploreShape.get_edges(section.shape)
        wires = WiresByConnectedEdges(edges).wires
        w = LengthOfShapes(wires).longest_shape
        w = CheckShape.to_wire(w)
        cref = ExploreShape.curve_of_shape(w)

        # Orient cref so that p1 is nearest (u1, v1) on the body
        umin, vmin = body.sref.u1, body.sref.v1
        p0 = body.eval(umin, vmin)
        p1 = cref.eval(cref.u1)
        p2 = cref.eval(cref.u2)
        d1 = p0.distance(p1)
        d2 = p0.distance(p2)
        if d2 < d1:
            cref.reverse()

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._spar = Spar(label, shape, cref, srf, group)

    @property
    def spar(self):
        """
        :return: The spar.
        :rtype: afem.structure.entities.Spar
        """
        return self._spar


class SparByShape(object):
    """
    Create a spar using a basis shape.

    :param str label: Part label.
    :param OCCT.TopoDS.TopoDS_Shape basis_shape: The basis shape.
    :param afem.oml.entities.Body body: The body.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        * The reference curve of the part will be untrimmed in this method. It
          will be determined by the intersection between *shape* and the body
          reference shape.

        * The reference surface of the part will be taken as the underlying
          surface of the largest face in the basis shape.
    """

    def __init__(self, label, basis_shape, body, group=None):
        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape,
                                  approximate=True)
        edges = ExploreShape.get_edges(section.shape)
        wires = WiresByConnectedEdges(edges).wires
        w = LengthOfShapes(wires).longest_shape
        w = CheckShape.to_wire(w)
        cref = ExploreShape.curve_of_shape(w)

        # Orient cref so that p1 is nearest (u1, v1) on the body
        umin, vmin = body.sref.u1, body.sref.v1
        p0 = body.eval(umin, vmin)
        p1 = cref.eval(cref.u1)
        p2 = cref.eval(cref.u2)
        d1 = p0.distance(p1)
        d2 = p0.distance(p2)
        if d2 < d1:
            cref.reverse()

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        basis_shape = common.shape

        # Get reference surface
        sref = ExploreShape.surface_of_shape(basis_shape)

        # Create the part
        self._spar = Spar(label, basis_shape, cref, sref, group)

    @property
    def spar(self):
        """
        :return: The spar.
        :rtype: afem.structure.entities.Spar
        """
        return self._spar


class SparBetweenShapes(SparByPoints):
    """
    Create a spar between shapes. This method uses the shapes to define
    starting and ending points by intersecting the reference curve and then
    uses :class:`.SparByPoints`.

    :param str label: Part label.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.body body: The body.
    :param basis_shape: The basis shape.
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, shape1, shape2, body, basis_shape, group=None):
        if isinstance(basis_shape, Surface):
            shape = FaceBySurface(basis_shape).face
        else:
            shape = basis_shape

        shape1 = shape_of_entity(shape1)
        shape2 = shape_of_entity(shape2)

        wing_basis_edges = IntersectShapes(shape, body.sref_shape).shape
        p1_shape = IntersectShapes(shape1, wing_basis_edges).shape
        p2_shape = IntersectShapes(shape2, wing_basis_edges).shape
        v1 = ExploreShape.get_vertices(p1_shape)[0]
        v2 = ExploreShape.get_vertices(p2_shape)[0]
        p1 = ExploreShape.pnt_of_vertex(v1)
        p2 = ExploreShape.pnt_of_vertex(v2)

        super(SparBetweenShapes, self).__init__(label, p1, p2, body,
                                                basis_shape, group)


class SparsBetweenPlanesByNumber(object):
    """
    Create a specified number of planar spars between two planes. This method
    uses :class:`.PlanesBetweenPlanesByNumber` and :class:`.SparBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, pln1, pln2, n, shape1, shape2, body, d1=None,
                 d2=None, first_index=1, delimiter=' ', group=None):
        n = int(n)
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByNumber(pln1, pln2, n, d1, d2)

        self._spars = []
        self._nspars = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            spar = SparBetweenShapes(label_indx, shape1, shape2, body,
                                     basis_shape, group).spar
            first_index += 1
            self._spars.append(spar)
        self._next_index = first_index

    @property
    def nspars(self):
        """
        :return: The number of spars.
        :rtype: int
        """
        return self._nspars

    @property
    def spars(self):
        """
        :return: The spars.
        :rtype: list[afem.structure.entities.Spar]
        """
        return self._spars

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class SparsBetweenPlanesByDistance(object):
    """
    Create planar spars between two planes using a maximum spacing. This method
    uses :class:`.PlanesBetweenPlanesByDistance` and
    :class:`.SparBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param float maxd: The maximum allowed spacing. The actual spacing will
        be adjusted to not to exceed this value.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of parts to create.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, pln1, pln2, maxd, shape1, shape2, body, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' ', group=None):
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._spars = []
        self._nspars = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            spar = SparBetweenShapes(label_indx, shape1, shape2, body,
                                     basis_shape, group).spar
            first_index += 1
            self._spars.append(spar)
        self._next_index = first_index

    @property
    def nspars(self):
        """
        :return: The number of spars.
        :rtype: int
        """
        return self._nspars

    @property
    def spars(self):
        """
        :return: The spars.
        :rtype: list[afem.structure.entities.Spar]
        """
        return self._spars

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class SparsAlongCurveByNumber(object):
    """
    Create a specified number of planar spars along a curve. This method
    uses :class:`.PlanesAlongCurveByNumber` and :class:`.SparBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param int n: The number of parts.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param afem.geometry.entities.Plane ref_pln: The normal of this plane
        will be used to define the normal of all planes along the curve. If
        no plane is provided, then the first derivative of the curve will
        define the plane normal.
    :param float u1: The parameter of the first plane (default=crv.u1).
    :param float u2: The parameter of the last plane (default=crv.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param float tol: Tolerance.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, crv, n, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, first_index=1,
                 delimiter=' ', tol=1.0e-7, group=None):
        n = int(n)
        first_index = int(first_index)

        builder = PlanesAlongCurveByNumber(crv, n, ref_pln, u1, u2, d1, d2,
                                           tol)

        self._spars = []
        self._nspars = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            spar = SparBetweenShapes(label_indx, shape1, shape2, body,
                                     basis_shape, group).spar
            first_index += 1
            self._spars.append(spar)
        self._next_index = first_index

    @property
    def nspars(self):
        """
        :return: The number of spars.
        :rtype: int
        """
        return self._nspars

    @property
    def spars(self):
        """
        :return: The spars.
        :rtype: list[afem.structure.entities.Spar]
        """
        return self._spars

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class SparsAlongCurveByDistance(object):
    """
    Create planar spars along a curve using a maximum spacing. This method
    uses :class:`.PlanesAlongCurveByDistance` and :class:`.SparBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param float maxd: The maximum allowed spacing between planes. The
        actual spacing will be adjusted to not to exceed this value.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param afem.geometry.entities.Plane ref_pln: The normal of this plane
        will be used to define the normal of all planes along the curve. If
        no plane is provided, then the first derivative of the curve will
        define the plane normal.
    :param float u1: The parameter of the first plane (default=crv.u1).
    :param float u2: The parameter of the last plane (default=crv.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of planes to create.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param float tol: Tolerance.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, crv, maxd, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, nmin=0, first_index=1,
                 delimiter=' ', tol=1.0e-7, group=None):
        first_index = int(first_index)

        builder = PlanesAlongCurveByDistance(crv, maxd, ref_pln, u1, u2, d1,
                                             d2, nmin, tol)

        self._spars = []
        self._nspars = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            spar = SparBetweenShapes(label_indx, shape1, shape2, body,
                                     basis_shape, group).spar
            first_index += 1
            self._spars.append(spar)
        self._next_index = first_index

    @property
    def nspars(self):
        """
        :return: The number of spars.
        :rtype: int
        """
        return self._nspars

    @property
    def spars(self):
        """
        :return: The spars.
        :rtype: list[afem.structure.entities.Spar]
        """
        return self._spars

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


# RIB -------------------------------------------------------------------------

class RibByParameters(object):
    """
    Create a rib between wing parameters.

    :param str label: Part label.
    :param float u1: Starting point u-parameter.
    :param float v1: Starting point v-parameter.
    :param float u2: Ending point u-parameter.
    :param float v2: Ending point v-parameter.
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        * If a basis shape is provided, then the starting and ending parameters
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, u1, v1, u2, v2, body, basis_shape=None,
                 group=None):
        # Determine reference surface and basis shape
        if basis_shape is None:
            pln = body.extract_plane(u1, v1, u2, v2)
            sref = pln
            basis_shape = FaceBySurface(pln).face
        elif isinstance(basis_shape, Surface):
            sref = basis_shape
            basis_shape = FaceBySurface(sref).face
        else:
            sref = ExploreShape.surface_of_shape(basis_shape)

        # Extract cref
        cref = body.extract_curve(u1, v1, u2, v2, basis_shape)

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._rib = Rib(label, shape, cref, sref, group)

    @property
    def rib(self):
        """
        :return: The rib.
        :rtype: afem.structure.entities.Rib
        """
        return self._rib


class RibByPoints(RibByParameters):
    """
    Create a rib between two points. This method inverts the starting and
    ending points and then uses :class:`.RibByParameters`.

    :param str label: Part label.
    :param point_like p1: Starting point.
    :param point_like p2: End point.
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        * The starting and ending points should be on or near the body
          reference surface since they are projected to find parameters.

        * If a basis shape is provided, then the starting and ending points
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, p1, p2, body, basis_shape=None, group=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        # Invert points
        u1, v1 = body.invert(p1)
        u2, v2 = body.invert(p2)

        # Use SparByParameters
        super(RibByPoints, self).__init__(label, u1, v1, u2, v2, body,
                                          basis_shape, group)


class RibBySurface(object):
    """
    Create a rib using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *srf* and the body
        reference shape.
    """

    def __init__(self, label, srf, body, group=None):
        basis_shape = FaceBySurface(srf).face

        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape,
                                  approximate=True)
        edges = ExploreShape.get_edges(section.shape)
        wires = WiresByConnectedEdges(edges).wires
        w = LengthOfShapes(wires).longest_shape
        w = CheckShape.to_wire(w)
        cref = ExploreShape.curve_of_shape(w)

        # Orient cref so that p1 is nearest (u1, v1) on the body
        umin, vmin = body.sref.u1, body.sref.v1
        p0 = body.eval(umin, vmin)
        p1 = cref.eval(cref.u1)
        p2 = cref.eval(cref.u2)
        d1 = p0.distance(p1)
        d2 = p0.distance(p2)
        if d2 < d1:
            cref.reverse()

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._rib = Rib(label, shape, cref, srf, group)

    @property
    def rib(self):
        """
        :return: The rib.
        :rtype: afem.structure.entities.Rib
        """
        return self._rib


class RibByShape(object):
    """
    Create a rib using a basis shape.

    :param str label: Part label.
    :param OCCT.TopoDS.TopoDS_Shape basis_shape: The basis shape.
    :param afem.oml.entities.Body body: The body.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    .. note::

        * The reference curve of the part will be untrimmed in this method. It
          will be determined by the intersection between *shape* and the body
          reference shape.

        * The reference surface of the part will be taken as the underlying
          surface of the largest face in the basis shape.
    """

    def __init__(self, label, basis_shape, body, group=None):
        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape,
                                  approximate=True)
        edges = ExploreShape.get_edges(section.shape)
        wires = WiresByConnectedEdges(edges).wires
        w = LengthOfShapes(wires).longest_shape
        w = CheckShape.to_wire(w)
        cref = ExploreShape.curve_of_shape(w)

        # Orient cref so that p1 is nearest (u1, v1) on the body
        umin, vmin = body.sref.u1, body.sref.v1
        p0 = body.eval(umin, vmin)
        p1 = cref.eval(cref.u1)
        p2 = cref.eval(cref.u2)
        d1 = p0.distance(p1)
        d2 = p0.distance(p2)
        if d2 < d1:
            cref.reverse()

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        basis_shape = common.shape

        # Get reference surface
        sref = ExploreShape.surface_of_shape(basis_shape)

        # Create the part
        self._rib = Rib(label, basis_shape, cref, sref, group)

    @property
    def rib(self):
        """
        :return: The rib.
        :rtype: afem.structure.entities.Rib
        """
        return self._rib


class RibBetweenShapes(RibByPoints):
    """
    Create a rib between shapes. This method uses the shapes to define
    starting and ending points by intersecting the reference curve and then
    uses :class:`.RibByPoints`.

    :param str label: Part label.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape.
    :type basis_shape: afem.geometry.entities.Surface or
        OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, shape1, shape2, body, basis_shape, group=None):
        if isinstance(basis_shape, Surface):
            shape = FaceBySurface(basis_shape).face
        else:
            shape = basis_shape

        shape1 = shape_of_entity(shape1)
        shape2 = shape_of_entity(shape2)

        wing_basis_edges = IntersectShapes(shape, body.sref_shape).shape
        p1_shape = IntersectShapes(shape1, wing_basis_edges).shape
        p2_shape = IntersectShapes(shape2, wing_basis_edges).shape
        v1 = ExploreShape.get_vertices(p1_shape)[0]
        v2 = ExploreShape.get_vertices(p2_shape)[0]
        p1 = ExploreShape.pnt_of_vertex(v1)
        p2 = ExploreShape.pnt_of_vertex(v2)

        super(RibBetweenShapes, self).__init__(label, p1, p2, body,
                                               basis_shape, group)


class RibByOrientation(RibBySurface):
    """
    Create a planar rib using rotation angles. This creates the plane then uses
    :class:`.RibBySurface`.

    :param str label: Part label.
    :param point_like origin: The origin of the plane after rotating and
        translating from the global origin.
    :param afem.oml.entities.Body body: The body.
    :param float alpha: Rotation in degrees about global x-axis.
    :param float beta: Rotation in degrees about global y-axis.
    :param float gamma: Rotation in degrees about global z-axis.
    :param str axes: The axes for the original plane before rotation and
        translation ('xy', 'xz', 'yz').
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, origin, body, alpha=0., beta=0., gamma=0.,
                 axes='xz', group=None):
        pln = PlaneByOrientation(origin, axes, alpha, beta, gamma).plane

        super(RibByOrientation, self).__init__(label, pln, body, group)


class RibsBetweenPlanesByNumber(object):
    """
    Create a specified number of planar ribs between two planes. This method
    uses :class:`.PlanesBetweenPlanesByNumber` and :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, pln1, pln2, n, shape1, shape2, body, d1=None,
                 d2=None, first_index=1, delimiter=' ', group=None):
        n = int(n)
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByNumber(pln1, pln2, n, d1, d2)

        self._ribs = []
        self._nribs = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            rib = RibBetweenShapes(label_indx, shape1, shape2, body,
                                   basis_shape, group).rib
            first_index += 1
            self._ribs.append(rib)
        self._next_index = first_index

    @property
    def nribs(self):
        """
        :return: The number of ribs.
        :rtype: int
        """
        return self._nribs

    @property
    def ribs(self):
        """
        :return: The ribs.
        :rtype: list[afem.structure.entities.Rib]
        """
        return self._ribs

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class RibsBetweenPlanesByDistance(object):
    """
    Create planar ribs between two planes using a maximum spacing. This method
    uses :class:`.PlanesBetweenPlanesByDistance` and
    :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param float maxd: The maximum allowed spacing. The actual spacing will
        be adjusted to not to exceed this value.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of parts to create.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, pln1, pln2, maxd, shape1, shape2, body, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' ', group=None):
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._ribs = []
        self._nribs = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            rib = RibBetweenShapes(label_indx, shape1, shape2, body,
                                   basis_shape, group).rib
            first_index += 1
            self._ribs.append(rib)
        self._next_index = first_index

    @property
    def nribs(self):
        """
        :return: The number of ribs.
        :rtype: int
        """
        return self._nribs

    @property
    def ribs(self):
        """
        :return: The ribs.
        :rtype: list[afem.structure.entities.Rib]
        """
        return self._ribs

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class RibsAlongCurveByNumber(object):
    """
    Create a specified number of planar ribs along a curve. This method
    uses :class:`.PlanesAlongCurveByNumber` and :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param int n: The number of parts.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param afem.geometry.entities.Plane ref_pln: The normal of this plane
        will be used to define the normal of all planes along the curve. If
        no plane is provided, then the first derivative of the curve will
        define the plane normal.
    :param float u1: The parameter of the first plane (default=crv.u1).
    :param float u2: The parameter of the last plane (default=crv.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param float tol: Tolerance.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, crv, n, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, first_index=1,
                 delimiter=' ', tol=1.0e-7, group=None):
        n = int(n)
        first_index = int(first_index)

        builder = PlanesAlongCurveByNumber(crv, n, ref_pln, u1, u2, d1, d2,
                                           tol)

        self._ribs = []
        self._nribs = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            rib = RibBetweenShapes(label_indx, shape1, shape2, body,
                                   basis_shape, group).rib
            first_index += 1
            self._ribs.append(rib)
        self._next_index = first_index

    @property
    def nribs(self):
        """
        :return: The number of ribs.
        :rtype: int
        """
        return self._nribs

    @property
    def ribs(self):
        """
        :return: The ribs.
        :rtype: list[afem.structure.entities.Rib]
        """
        return self._ribs

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class RibsAlongCurveByDistance(object):
    """
    Create planar ribs along a curve using a maximum spacing. This method
    uses :class:`.PlanesAlongCurveByDistance` and :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param float maxd: The maximum allowed spacing between planes. The
        actual spacing will be adjusted to not to exceed this value.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param afem.geometry.entities.Plane ref_pln: The normal of this plane
        will be used to define the normal of all planes along the curve. If
        no plane is provided, then the first derivative of the curve will
        define the plane normal.
    :param float u1: The parameter of the first plane (default=crv.u1).
    :param float u2: The parameter of the last plane (default=crv.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of planes to create.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param float tol: Tolerance.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, crv, maxd, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, nmin=0, first_index=1,
                 delimiter=' ', tol=1.0e-7, group=None):
        first_index = int(first_index)

        builder = PlanesAlongCurveByDistance(crv, maxd, ref_pln, u1, u2, d1,
                                             d2, nmin, tol)

        self._ribs = []
        self._nribs = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            rib = RibBetweenShapes(label_indx, shape1, shape2, body,
                                   basis_shape, group).rib
            first_index += 1
            self._ribs.append(rib)
        self._next_index = first_index

    @property
    def nribs(self):
        """
        :return: The number of ribs.
        :rtype: int
        """
        return self._nribs

    @property
    def ribs(self):
        """
        :return: The ribs.
        :rtype: list[afem.structure.entities.Rib]
        """
        return self._ribs

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class RibsAlongCurveAndSurfaceByDistance(object):
    """
    Create planar ribs along a curve and surface using a maximum spacing. This
    method uses :class:`.PlanesAlongCurveAndSurfaceByDistance` and
    :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param afem.geometry.entities.Surface srf: The surface.
    :param float maxd: The maximum allowed spacing between planes. The
        actual spacing will be adjusted to not to exceed this value.
    :param shape1: Starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param shape2: Ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface or afem.structure.entities.Part
    :param afem.oml.entities.Body body: The body.
    :param float u1: The parameter of the first plane (default=crv.u1).
    :param float u2: The parameter of the last plane (default=crv.u2).
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param float rot_x: The rotation angles of each plane along their local
        x-axis in degrees.
    :param float rot_y: The rotation angles of each plane along their local
        y-axis in degrees.
    :param int nmin: Minimum number of planes to create.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param float tol: Tolerance.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, crv, srf, maxd, shape1, shape2, body,
                 u1=None, u2=None, d1=None, d2=None, rot_x=None, rot_y=None,
                 nmin=0, first_index=1, delimiter=' ', tol=1.0e-7, group=None):
        first_index = int(first_index)

        builder = PlanesAlongCurveAndSurfaceByDistance(crv, srf, maxd, u1, u2,
                                                       d1, d2, nmin, tol)
        if rot_x is not None:
            builder.rotate_x(rot_x)
        if rot_y is not None:
            builder.rotate_y(rot_y)

        self._ribs = []
        self._nribs = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            rib = RibBetweenShapes(label_indx, shape1, shape2, body,
                                   basis_shape, group).rib
            first_index += 1
            self._ribs.append(rib)
        self._next_index = first_index

    @property
    def nribs(self):
        """
        :return: The number of ribs.
        :rtype: int
        """
        return self._nribs

    @property
    def ribs(self):
        """
        :return: The ribs.
        :rtype: list[afem.structure.entities.Rib]
        """
        return self._ribs

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


# BULKHEAD --------------------------------------------------------------------

class BulkheadBySurface(object):
    """
    Create a bulkhead using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, srf, body, group=None):
        basis_shape = FaceBySurface(srf).face

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._bh = Bulkhead(label, shape, None, srf, group)

    @property
    def bulkhead(self):
        """
        :return: The bulkhead.
        :rtype: afem.structure.entities.Bulkhead
        """
        return self._bh


# FLOOR -----------------------------------------------------------------------

class FloorBySurface(object):
    """
    Create a floor using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, srf, body, group=None):
        basis_shape = FaceBySurface(srf).face

        # Build part shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._floor = Floor(label, shape, None, srf, group)

    @property
    def floor(self):
        """
        :return: The floor.
        :rtype: afem.structure.entities.Floor
        """
        return self._floor


# FRAME -----------------------------------------------------------------------

class FrameByPlane(object):
    """
    Create a frame using a plane. A plane is required since the shape of
    frame is formed using an offset algorithm that requires planar shapes.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln: The plane.
    :param afem.oml.entities.Body body: The body.
    :param float height: The height. The absolute value is used.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, pln, body, height, group=None):
        basis_shape = FaceBySurface(pln).face

        # Find initial shape
        common = CommonShapes(basis_shape, body.solid)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Get outer (free) edge of shape which should be a closed wire. Use
        #  the longest wire if necessary.
        closed_wires = ExploreFreeEdges(shape).closed_wires
        if len(closed_wires) > 1:
            outer_wire = LengthOfShapes(closed_wires).longest_shape
        else:
            outer_wire = closed_wires[0]

        # Offset the outer wire and concatenate it
        inner_wire = WireByPlanarOffset(outer_wire, -abs(height)).wire
        inner_wire = WireByConcat(inner_wire).wire

        # Create inner and outer shapes and cut one from the other
        inner_face = FaceByPlanarWire(inner_wire).face
        cut = CutShapes(shape, inner_face)
        if not cut.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = cut.shape

        # Create the part
        self._frame = Frame(label, shape, None, pln, group)

    @property
    def frame(self):
        """
        :return: The frame.
        :rtype: afem.structure.entities.Frame
        """
        return self._frame


class FramesByPlanes(object):
    """
    Create frames using a list of planes. This method uses
    :class:`.FrameByPlane`.

    :param str label: Part label.
    :param list[afem.geometry.entities.Plane] plns: The planes.
    :param afem.oml.entities.Body body: The body.
    :param float height: The height.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, plns, body, height, first_index=1,
                 delimiter=' ', group=None):
        first_index = int(first_index)

        self._frames = []
        self._nframes = len(plns)
        for pln in plns:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, body, height, group).frame
            first_index += 1
            self._frames.append(frame)
        self._next_index = first_index

    @property
    def nframes(self):
        """
        :return: The number of frames.
        :rtype: int
        """
        return self._nframes

    @property
    def frames(self):
        """
        :return: The frames.
        :rtype: list[afem.structure.entities.Frame]
        """
        return self._frames

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class FramesBetweenPlanesByNumber(object):
    """
    Create a specified number of frames between two planes. This method uses
    :class:`.PlanesBetweenPlanesByNumber` and :class:`.FrameByPlane`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param afem.oml.entities.Body body: The body.
    :param float height: The height.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, pln1, pln2, n, body, height, d1=None,
                 d2=None, first_index=1, delimiter=' ', group=None):
        n = int(n)
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByNumber(pln1, pln2, n, d1, d2)

        self._frames = []
        self._nframes = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, body, height, group).frame
            first_index += 1
            self._frames.append(frame)
        self._next_index = first_index

    @property
    def nframes(self):
        """
        :return: The number of frames.
        :rtype: int
        """
        return self._nframes

    @property
    def frames(self):
        """
        :return: The frames.
        :rtype: list[afem.structure.entities.Frame]
        """
        return self._frames

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


class FramesBetweenPlanesByDistance(object):
    """
    Create frames between two planes using a maximum spacing. This method uses
    :class:`.PlanesBetweenPlanesByDistance` and :class:`.FrameByPlane`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param float maxd: The maximum allowed spacing. The actual spacing will
        be adjusted to not to exceed this value.
    :param afem.oml.entities.Body body: The body.
    :param float height: The height.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int nmin: Minimum number of parts to create.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, pln1, pln2, maxd, body, height, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' ', group=None):
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._frames = []
        self._nframes = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, body, height, group).frame
            first_index += 1
            self._frames.append(frame)
        self._next_index = first_index

    @property
    def nframes(self):
        """
        :return: The number of frames.
        :rtype: int
        """
        return self._nframes

    @property
    def frames(self):
        """
        :return: The frames.
        :rtype: list[afem.structure.entities.Frame]
        """
        return self._frames

    @property
    def spacing(self):
        """
        :return: The spacing between first and second parts.
        :rtype: float
        """
        return self._ds

    @property
    def next_index(self):
        """
        :return: The next index.
        :rtype: int
        """
        return self._next_index


# SKIN ------------------------------------------------------------------------

class SkinBySolid(object):
    """
    Create a skin part from the outer shell of the solid.

    :param str label: Part label.
    :param OCCT.TopoDS.TopoDS_Solid solid: The solid.
    :param bool copy: Option to copy the outer shell.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, solid, copy=False, group=None):
        shell = ExploreShape.outer_shell(solid)
        if copy:
            shell = ExploreShape.copy_shape(shell, False)

        self._skin = Skin(label, shell, None, None, group)

    @property
    def skin(self):
        """
        :return: The skin.
        :rtype: afem.structure.entities.Skin
        """
        return self._skin


class SkinByBody(SkinBySolid):
    """
    Create a skin part from the outer shell of a body.

    :param str label: Part label.
    :param afem.oml.entities.Body body: The body.
    :param bool copy: Option to copy the outer shell.
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, body, copy=False, group=None):
        super(SkinByBody, self).__init__(label, body.solid, copy, group)


# STRINGER --------------------------------------------------------------------

class StringerByShape(object):
    """
    Create a Stringer using a basis shape.

    :param str label: Part label.
    :param basis_shape: The basis shape that will define the path of the
        stringer.
    :type basis_shape: OCCT.TopoDS.TopoDS_Shape or
        afem.geometry.entities.Surface
    :param support_shape: The shape that will defines the normal direction along
        the path.
    :type support_shape: OCCT.TopoDS.TopoDS_Shape or
        afem.geometry.entities.Surface
    :param float height: The height.
    :param float angle: The runout angle at each end.
    :param shape1: The starting shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param shape2: The ending shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, basis_shape, support_shape, height, angle=30.,
                 shape1=None, shape2=None, group=None):

        # Convert to shapes
        support_shape = CheckShape.to_shape(support_shape)
        shape1 = CheckShape.to_shape(shape1)
        shape2 = CheckShape.to_shape(shape2)

        # Intersect to find spine
        spine = IntersectShapes(support_shape, basis_shape, True).shape

        # Build the wire
        tool = WiresByShape(spine)
        if tool.nwires == 0:
            msg = 'Failed to find a wire for the stringer.'
            raise RuntimeError(msg)
        elif tool.nwires > 1:
            msg = 'Found more than one wire for the stringer. Using longest.'
            warn(msg, RuntimeWarning)
            spine = LengthOfShapes(tool.wires).longest_shape
        else:
            spine = tool.wires[0]

        # Trim spine between shapes
        if (shape1 and shape2) is not None:
            spine = TrimOpenWire(spine, shape1, shape2).trimmed_wire
        elif shape1 is not None and shape2 is None:
            spine = TrimOpenWire(spine, shape1).last_wire
        elif shape1 is None and shape2 is not None:
            spine = TrimOpenWire(spine, shape2).first_wire

        # Determine length along spine given height and angle
        dx = height / tan(radians(angle))

        # Make the profiles
        adp_crv = BRepAdaptor_CompCurve(spine, True)

        p0 = CheckGeom.to_point(adp_crv.Value(dx))
        dss = DistanceShapeToShape(support_shape, p0)
        vn = CheckGeom.to_vector(dss.normal_on_shape1(1))
        vn.scale(height)
        p1 = p0.copy()
        p1.translate(vn)
        profile1 = EdgeByPoints(p0, p1).edge

        p0 = CheckGeom.to_point(adp_crv.Value(adp_crv.LastParameter() - dx))
        dss = DistanceShapeToShape(support_shape, p0)
        vn = CheckGeom.to_vector(dss.normal_on_shape1(1))
        vn.scale(height)
        p1 = p0.copy()
        p1.translate(vn)
        profile2 = EdgeByPoints(p0, p1).edge

        # Trim spine at the first and last profiles
        trim = TrimOpenWire(spine, profile1, profile2)
        first_wire = trim.first_wire
        middle_wire = trim.trimmed_wire
        last_wire = trim.last_wire

        # Make first segment
        pipe = SweepShapeWithNormal(first_wire, support_shape)
        pipe.add_profile(trim.all_vertices[0], True)
        pipe.add_profile(profile1, True)
        pipe.build()
        first_shape = pipe.shape

        # Make middle segment
        pipe = SweepShapeWithNormal(middle_wire, support_shape)
        pipe.add_profile(profile1, True)
        pipe.add_profile(profile2, True)
        pipe.build()
        middle_shape = pipe.shape

        # Make last segment
        pipe = SweepShapeWithNormal(last_wire, support_shape)
        pipe.add_profile(profile2, True)
        pipe.add_profile(trim.all_vertices[-1], True)
        pipe.build()
        last_shape = pipe.shape

        # Sew segments
        cmp = CompoundByShapes(
            [first_shape, middle_shape, last_shape]).compound
        shape = SewShape(cmp).sewed_shape

        self._stringer = Stringer(label, shape, group=group)

    @property
    def stringer(self):
        """
        :return: The stringer.
        :rtype: afem.structure.entities.Stringer
        """
        return self._stringer


# BEAM2D ----------------------------------------------------------------------

class Beam2DBySweep(object):
    """
    Create a Beam2D by sweeping a profile along a path.

    :param str label: Part label.
    :param spine: The path for the sweep.
    :type spine: afem.geometry.entities.Curve or OCCT.TopoDS.TopoDS_Edge or
        OCCT.TopoDS.TopoDS_Wire
    :param profile: The profile to sweep.
    :type profile: afem.geometry.entities.Curve or OCCT.TopoDS.TopoDS_Shape
    :param group: The group to add the part to. If not provided the part will
        be added to the active group.
    :type group: str or afem.structure.group.Group or None
    """

    def __init__(self, label, spine, profile, group=None):
        cref = None
        if CheckGeom.is_curve(spine):
            cref = spine
        elif CheckShape.is_edge(spine):
            cref = ExploreShape.curve_of_edge(spine)

        spine = CheckShape.to_wire(spine)
        profile = CheckShape.to_shape(profile)

        if cref is None:
            cref = ExploreShape.curve_of_wire(spine)

        tool = SweepShape(spine, profile)
        self._beam = Beam2D(label, tool.shape, cref, None, group)

    @property
    def beam2d(self):
        """
        :return: The beam.
        :rtype: afem.structure.entities.Beam2D
        """
        return self._beam


if __name__ == "__main__":
    import doctest

    doctest.testmod()
