from OCC.TopoDS import TopoDS_Shape, TopoDS_Solid

from afem.geometry.check import CheckGeom
from afem.geometry.create import *
from afem.geometry.entities import *
from afem.oml.entities import Body
from afem.structure.entities import *
from afem.topology.bop import *
from afem.topology.check import CheckShape
from afem.topology.create import *
from afem.topology.explore import ExploreFreeEdges, ExploreShape
from afem.topology.props import *

__all__ = ["CurvePartByShape", "BeamByShape", "BeamByCurve", "BeamByPoints",
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
           "SkinByBody"]


def _validate_type(input_, type_, required=True):
    """
    Validate the input is of the required type.
    """
    if not required and not input_:
        return None
    if not isinstance(input_, type_):
        msg = 'Invalid input type.'
        raise TypeError(msg)


# CURVE PART ------------------------------------------------------------------

class CurvePartByShape(object):
    """
    Create a curve part using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: afem.geometry.entities.Curve or OCC.TopoDS.TopoDS_Shape
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> from afem.structure import CurvePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> part = CurvePartByShape('part', e).curve_part
    """

    def __init__(self, label, shape, cref=None, assy=None):
        _validate_type(shape, (Curve, TopoDS_Shape))
        _validate_type(cref, Curve, False)

        if isinstance(shape, Curve):
            cref = shape
            shape = EdgeByCurve(shape).edge

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._curve_part = CurvePart(label, shape, cref, assy)

    @property
    def curve_part(self):
        """
        :return: The curve part.
        :rtype: afem.structure.entities.CurvePart
        """
        return self._curve_part


# BEAM -----------------------------------------------------------------------

class BeamByShape(object):
    """
    Create a beam using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: afem.geometry.entities.Curve or OCC.TopoDS.TopoDS_Shape
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> from afem.structure import BeamByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> beam = BeamByShape('part', e).beam
    """

    def __init__(self, label, shape, cref=None, assy=None):
        _validate_type(shape, (Curve, TopoDS_Shape))
        _validate_type(cref, Curve, False)

        if isinstance(shape, Curve):
            shape = EdgeByCurve(shape).edge

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._beam = Beam(label, shape, cref, assy)

    @property
    def beam(self):
        """
        :return: The beam.
        :rtype: afem.structure.entities.Beam
        """
        return self._beam


class BeamByCurve(BeamByShape):
    """
    Create a beam using a curve.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
    >>> from afem.structure import BeamByCurve
    >>> c = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> part = BeamByCurve('part', c).beam
    """

    def __init__(self, label, crv, assy=None):
        _validate_type(crv, Curve)

        e = EdgeByCurve(crv).edge
        super(BeamByCurve, self).__init__(label, e, crv, assy)


class BeamByPoints(BeamByCurve):
    """
    Create a beam between two points.

    :param str label: Part label.
    :param point_like p1: First point.
    :param point_like p2: Second point.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    Usage:

    >>> from afem.structure import BeamByPoints
    >>> beam = BeamByPoints('part', (0., 0., 0.), (10., 0., 0.)).beam
    >>> beam.length
    10.0
    """

    def __init__(self, label, p1, p2, assy=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        _validate_type(p1, Point)
        _validate_type(p2, Point)

        c = NurbsCurveByPoints([p1, p2]).curve
        super(BeamByPoints, self).__init__(label, c, assy)


# SURFACE PART ----------------------------------------------------------------

class SurfacePartByShape(object):
    """
    Create a surface part using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: afem.geometry.entities.Surface or OCC.TopoDS.TopoDS_Shape
    :param afem.geometry.entities.Curve cref: The reference curve.
    :param afem.geometry.entities.Surface sref: The reference surface. If
        not provided then the basis surface of the largest face of the shape
        will be used.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    Usage:

    >>> from afem.topology import EdgeByPoints, FaceByDrag
    >>> from afem.structure import SurfacePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> f = FaceByDrag(e, (0., 10., 0.)).face
    >>> part = SurfacePartByShape('part', f).surface_part
    >>> part.area
    99.99999999999997
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        _validate_type(shape, (Surface, TopoDS_Shape))
        _validate_type(cref, Curve, False)
        _validate_type(sref, Surface, False)

        if isinstance(shape, Surface):
            sref = shape
            shape = FaceBySurface(shape).face

        if sref is None:
            sref = ExploreShape.surface_of_shape(shape)

        self._surface_part = SurfacePart(label, shape, cref, sref, assy)

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
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.

    .. note::

        * If a basis shape is provided, then the starting and ending parameters
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, u1, v1, u2, v2, body, basis_shape=None,
                 assy=None):
        _validate_type(body, Body)
        _validate_type(basis_shape, (Surface, TopoDS_Shape), False)

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
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._spar = Spar(label, shape, cref, sref, assy)

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
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        * The starting and ending points should be on or near the body
          reference surface since they are projected to find parameters.

        * If a basis shape is provided, then the starting and ending points
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, p1, p2, body, basis_shape=None, assy=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        _validate_type(p1, Point)
        _validate_type(p2, Point)

        # Invert points
        u1, v1 = body.invert(p1)
        u2, v2 = body.invert(p2)

        # Use SparByParameters
        super(SparByPoints, self).__init__(label, u1, v1, u2, v2, body,
                                           basis_shape, assy)


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
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    :raise TypeError: If *e1* or *e2* are not *point_like* or a sequence of
        two surface parameters.
    """

    def __init__(self, label, e1, e2, body, basis_shape=None, assy=None):
        if len(e1) == 2:
            u1, v1 = e1
        elif CheckGeom.is_point_like(e1):
            u1, v1 = body.invert(e1)
        else:
            raise TypeError('Invalid type for first e1')

        if len(e2) == 2:
            u2, v2 = e2
        elif CheckGeom.is_point_like(e2):
            u2, v2 = body.invert(e2)
        else:
            raise TypeError('Invalid type for e2.')

        super(SparByEnds, self).__init__(label, u1, v1, u2, v2, body,
                                         basis_shape, assy)


class SparBySurface(object):
    """
    Create a spar using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *srf* and the body
        reference shape.
    """

    def __init__(self, label, srf, body, assy=None):
        _validate_type(srf, Surface)
        _validate_type(body, Body)

        basis_shape = FaceBySurface(srf).face

        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape)
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
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._spar = Spar(label, shape, cref, srf, assy)

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
    :param OCC.TopoDS.TopoDS_Shape basis_shape: The basis shape.
    :param afem.oml.entities.Body body: The body.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        * The reference curve of the part will be untrimmed in this method. It
          will be determined by the intersection between *shape* and the body
          reference shape.

        * The reference surface of the part will be taken as the underlying
          surface of the largest face in the basis shape.
    """

    def __init__(self, label, basis_shape, body, assy=None):
        _validate_type(basis_shape, TopoDS_Shape)
        _validate_type(body, Body)

        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape)
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
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        basis_shape = common.shape

        # Get reference surface
        sref = ExploreShape.surface_of_shape(basis_shape)

        # Create the part
        self._spar = Spar(label, basis_shape, cref, sref, assy)

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
    :type shape1: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface
    :param shape2: Ending shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface
    :param afem.oml.entities.body body: The body.
    :param basis_shape: The basis shape.
    :type basis_shape: afem.geometry.entities.Surface or
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, shape1, shape2, body, basis_shape, assy=None):
        _validate_type(shape1, (TopoDS_Shape, Curve, Surface))
        _validate_type(shape2, (TopoDS_Shape, Curve, Surface))
        _validate_type(body, Body)
        _validate_type(basis_shape, (Surface, TopoDS_Shape), False)

        if isinstance(basis_shape, Surface):
            shape = FaceBySurface(basis_shape).face
        else:
            shape = basis_shape

        shape1 = CheckShape.to_shape(shape1)
        shape2 = CheckShape.to_shape(shape2)

        wing_basis_edges = IntersectShapes(shape, body.sref_shape).shape
        p1_shape = IntersectShapes(shape1, wing_basis_edges).shape
        p2_shape = IntersectShapes(shape2, wing_basis_edges).shape
        v1 = ExploreShape.get_vertices(p1_shape)[0]
        v2 = ExploreShape.get_vertices(p2_shape)[0]
        p1 = ExploreShape.pnt_of_vertex(v1)
        p2 = ExploreShape.pnt_of_vertex(v2)

        super(SparBetweenShapes, self).__init__(label, p1, p2, body,
                                                basis_shape, assy)


class SparsBetweenPlanesByNumber(object):
    """
    Create a specified number of planar spars between two planes. This method
    uses :class:`.PlanesBetweenPlanesByNumber` and :class:`.SparBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, pln1, pln2, n, shape1, shape2, body, d1=None,
                 d2=None, first_index=1, delimiter=' ', assy=None):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(body, Body)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)

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
                                     basis_shape, assy).spar
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, pln1, pln2, maxd, shape1, shape2, body, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' ', assy=None):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(body, Body)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)

        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._spars = []
        self._nspars = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            spar = SparBetweenShapes(label_indx, shape1, shape2, body,
                                     basis_shape, assy).spar
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, crv, n, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, first_index=1,
                 delimiter=' ', tol=1.0e-7, assy=None):
        _validate_type(crv, Curve)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(body, Body)
        _validate_type(ref_pln, Plane, False)

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
                                     basis_shape, assy).spar
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, crv, maxd, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, nmin=0, first_index=1,
                 delimiter=' ', tol=1.0e-7, assy=None):
        _validate_type(crv, Curve)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(body, Body)
        _validate_type(ref_pln, Plane, False)

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
                                     basis_shape, assy).spar
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
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        * If a basis shape is provided, then the starting and ending parameters
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, u1, v1, u2, v2, body, basis_shape=None,
                 assy=None):
        _validate_type(body, Body)
        _validate_type(basis_shape, (Surface, TopoDS_Shape), False)

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
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._rib = Rib(label, shape, cref, sref, assy)

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
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        * The starting and ending points should be on or near the body
          reference surface since they are projected to find parameters.

        * If a basis shape is provided, then the starting and ending points
          should be near the intersection between this shape and the body
          reference shape.
    """

    def __init__(self, label, p1, p2, body, basis_shape=None, assy=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        _validate_type(p1, Point)
        _validate_type(p2, Point)

        # Invert points
        u1, v1 = body.invert(p1)
        u2, v2 = body.invert(p2)

        # Use SparByParameters
        super(RibByPoints, self).__init__(label, u1, v1, u2, v2, body,
                                          basis_shape, assy)


class RibBySurface(object):
    """
    Create a rib using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *srf* and the body
        reference shape.
    """

    def __init__(self, label, srf, body, assy=None):
        _validate_type(srf, Surface)
        _validate_type(body, Body)

        basis_shape = FaceBySurface(srf).face

        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape)
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
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._rib = Rib(label, shape, cref, srf, assy)

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
    :param OCC.TopoDS.TopoDS_Shape basis_shape: The basis shape.
    :param afem.oml.entities.Body body: The body.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    .. note::

        * The reference curve of the part will be untrimmed in this method. It
          will be determined by the intersection between *shape* and the body
          reference shape.

        * The reference surface of the part will be taken as the underlying
          surface of the largest face in the basis shape.
    """

    def __init__(self, label, basis_shape, body, assy=None):
        _validate_type(basis_shape, TopoDS_Shape)
        _validate_type(body, Body)

        # Build reference curve
        section = IntersectShapes(basis_shape, body.sref_shape)
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
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        basis_shape = common.shape

        # Get reference surface
        sref = ExploreShape.surface_of_shape(basis_shape)

        # Create the part
        self._rib = Rib(label, basis_shape, cref, sref, assy)

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
    :type shape1: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface
    :param shape2: Ending shape.
    :type shape2: OCC.TopoDS.TopoDS_Shape or afem.geometry.entities.Curve or
        afem.geometry.entities.Surface
    :param afem.oml.entities.Body body: The body.
    :param basis_shape: The basis shape.
    :type basis_shape: afem.geometry.entities.Surface or
        OCC.TopoDS.TopoDS_Shape
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, shape1, shape2, body, basis_shape, assy=None):
        _validate_type(shape1, (TopoDS_Shape, Curve, Surface))
        _validate_type(shape2, (TopoDS_Shape, Curve, Surface))
        _validate_type(body, Body)
        _validate_type(basis_shape, (Surface, TopoDS_Shape), False)

        if isinstance(basis_shape, Surface):
            shape = FaceBySurface(basis_shape).face
        else:
            shape = basis_shape

        shape1 = CheckShape.to_shape(shape1)
        shape2 = CheckShape.to_shape(shape2)

        wing_basis_edges = IntersectShapes(shape, body.sref_shape).shape
        p1_shape = IntersectShapes(shape1, wing_basis_edges).shape
        p2_shape = IntersectShapes(shape2, wing_basis_edges).shape
        v1 = ExploreShape.get_vertices(p1_shape)[0]
        v2 = ExploreShape.get_vertices(p2_shape)[0]
        p1 = ExploreShape.pnt_of_vertex(v1)
        p2 = ExploreShape.pnt_of_vertex(v2)

        super(RibBetweenShapes, self).__init__(label, p1, p2, body,
                                               basis_shape, assy)


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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, origin, body, alpha=0., beta=0., gamma=0.,
                 axes='xz', assy=None):
        pln = PlaneByOrientation(origin, axes, alpha, beta, gamma).plane

        super(RibByOrientation, self).__init__(label, pln, body, assy)


class RibsBetweenPlanesByNumber(object):
    """
    Create a specified number of planar ribs between two planes. This method
    uses :class:`.PlanesBetweenPlanesByNumber` and :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.entities.Plane pln1: The first plane.
    :param afem.geometry.entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, pln1, pln2, n, shape1, shape2, body, d1=None,
                 d2=None, first_index=1, delimiter=' ', assy=None):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(body, Body)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)

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
                                   basis_shape, assy).rib
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, pln1, pln2, maxd, shape1, shape2, body, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' ', assy=None):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(body, Body)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)

        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._ribs = []
        self._nribs = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            basis_shape = FaceBySurface(pln).face
            label_indx = delimiter.join([label, str(first_index)])
            rib = RibBetweenShapes(label_indx, shape1, shape2, body,
                                   basis_shape, assy).rib
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, crv, n, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, first_index=1,
                 delimiter=' ', tol=1.0e-7, assy=None):
        _validate_type(crv, Curve)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(body, Body)
        _validate_type(ref_pln, Plane, False)

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
                                   basis_shape, assy).rib
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, crv, maxd, shape1, shape2, body, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, nmin=0, first_index=1,
                 delimiter=' ', tol=1.0e-7, assy=None):
        _validate_type(crv, Curve)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(body, Body)
        _validate_type(ref_pln, Plane, False)

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
                                   basis_shape, assy).rib
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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, crv, srf, maxd, shape1, shape2, body,
                 u1=None, u2=None, d1=None, d2=None, rot_x=None, rot_y=None,
                 nmin=0, first_index=1, delimiter=' ', tol=1.0e-7, assy=None):
        _validate_type(crv, Curve)
        _validate_type(srf, Surface)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(body, Body)

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
                                   basis_shape, assy).rib
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


# BULKHEAD -------------------------------------------------------------------

class BulkheadBySurface(object):
    """
    Create a bulkhead using a surface.

    :param str label: Part label.
    :param afem.geometry.entities.Surface srf: The surface.
    :param afem.oml.entities.Body body: The body.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, srf, body, assy=None):
        _validate_type(srf, Surface)
        _validate_type(body, Body)

        basis_shape = FaceBySurface(srf).face

        # Build part shape
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._bh = Bulkhead(label, shape, None, srf, assy)

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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, srf, body, assy=None):
        _validate_type(srf, Surface)
        _validate_type(body, Body)

        basis_shape = FaceBySurface(srf).face

        # Build part shape
        common = CommonShapes(basis_shape, body)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape

        # Create the part
        self._floor = Floor(label, shape, None, srf, assy)

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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, pln, body, height, assy=None):
        _validate_type(pln, Plane)
        _validate_type(body, Body)

        basis_shape = FaceBySurface(pln).face

        # Find initial shape
        common = CommonShapes(basis_shape, body)
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
        self._frame = Frame(label, shape, None, pln, assy)

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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, plns, body, height, first_index=1,
                 delimiter=' ', assy=None):
        for pln in plns:
            _validate_type(pln, Plane)
        _validate_type(body, Body)

        first_index = int(first_index)

        self._frames = []
        self._nframes = len(plns)
        for pln in plns:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, body, height, assy).frame
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, pln1, pln2, n, body, height, d1=None,
                 d2=None, first_index=1, delimiter=' ', assy=None):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(body, Body)

        n = int(n)
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByNumber(pln1, pln2, n, d1, d2)

        self._frames = []
        self._nframes = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, body, height, assy).frame
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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, pln1, pln2, maxd, body, height, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' ', assy=None):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(body, Body)

        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._frames = []
        self._nframes = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, body, height, assy).frame
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
    :param OCC.TopoDS.TopoDS_Solid solid: The solid.
    :param bool copy: Option to copy the outer shell.
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, solid, copy=False, assy=None):
        _validate_type(solid, TopoDS_Solid)

        shell = ExploreShape.outer_shell(solid)
        if copy:
            shell = ExploreShape.copy_shape(shell, False)

        self._skin = Skin(label, shell, None, None, assy)

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
    :param assy: The assembly to add the part to. If not provided the part will
        be added to the active assembly.
    :type assy: str or afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, body, copy=False, assy=None):
        _validate_type(body, Body)
        super(SkinByBody, self).__init__(label, body, copy, assy)


# STRINGER --------------------------------------------------------------------

class StringerBySpine(object):
    pass


class StringerBySection(object):
    pass


class StringersBySections(object):
    pass


if __name__ == "__main__":
    import doctest

    doctest.testmod()


# OLD METHODS FOR REFERENCE ---------------------------------------------------

# def create_stiffener1d(surface_part, stiffener, label):
#     """
#     Create a 1-D stiffener on a surface part.
#     """
#     if not isinstance(surface_part, SurfacePart):
#         return None
#
#     cref = None
#     if CheckGeom.is_curve_like(stiffener):
#         cref = stiffener
#
#     # If stiffener is already a Stiffener1D part, split the two parts.
#     # Otherwise, find the intersection and create the stiffener.
#     if not isinstance(stiffener, Stiffener1D):
#         shape = ShapeTools.to_shape(stiffener)
#         if not shape:
#             return None
#         edges = ShapeTools.bsection(surface_part, shape, 'edge')
#         if not edges:
#             return None
#         wires = ShapeTools.connect_edges(edges)
#         if not wires:
#             return None
#         if len(wires) == 1:
#             curve_shape = wires[0]
#             w = curve_shape
#         else:
#             curve_shape = ShapeTools.make_compound(wires)
#             w = ShapeTools.longest_wire(wires)
#
#         # Build reference curve from longest wire.
#         if not cref:
#             cref = ShapeTools.curve_of_wire(w)
#
#         # Build stiffener
#         stiffener = Stiffener1D(label, curve_shape, cref)
#
#     # Split the parts.
#     surface_part.split(stiffener)
#
#     # Add to surface part.
#     surface_part.add_subpart(stiffener.label, stiffener)
#
#     return stiffener
#
#
# def create_stiffener2d_by_wire(etype, label, surface_part, wire, h,
#                                runout_angle, sref=None, assy=None):
#     """
#     Create a 2-D stiffener on a surface part by a wire.
#     """
#     cref = None
#     if CheckGeom.is_curve_like(wire):
#         cref = wire
#
#     # Create spine.
#     spine0 = ShapeTools.to_wire(wire)
#
#     # Build reference curve from spine.
#     if not cref:
#         cref = ShapeTools.curve_of_wire(spine0)
#
#     # Calculate dx along spine to get proper run-out angle
#     dx = h / tan(radians(runout_angle))
#
#     # Find points at dx from each end.
#     adp_crv = BRepAdaptor_CompCurve(spine0)
#     abs_pnt = GCPnts_AbscissaPoint(adp_crv, dx, adp_crv.FirstParameter())
#     u = abs_pnt.Parameter()
#     p1 = adp_crv.Value(u)
#     abs_pnt = GCPnts_AbscissaPoint(adp_crv, -dx, adp_crv.LastParameter())
#     u = abs_pnt.Parameter()
#     p2 = adp_crv.Value(u)
#
#     # Create profiles normal to surface part.
#     profile1 = ShapeTools.shape_normal_profile(p1, surface_part, h)
#     profile2 = ShapeTools.shape_normal_profile(p2, surface_part, h)
#     if not profile1 or not profile2:
#         return None
#
#     # Create new spine between the profiles.
#     splitter = ShapeTools.make_compound([p1, p2])
#     spine1 = ShapeTools.split_wire(spine0, splitter)
#     if not spine1:
#         return None
#     # Create new spine omitting first and last edges.
#     ordered_edges = spine1.ordered_edges
#     if len(ordered_edges) < 3:
#         return None
#     try:
#         spine1 = ShapeTools.connect_edges(ordered_edges[1:-1])[0]
#     except IndexError:
#         return None
#
#     # Create main shape.
#     shape1 = ShapeTools.make_pipe_shell(spine1, [profile1, profile2],
#                                         surface_part, True)
#     if not shape1:
#         return None
#
#     # TODO Handle case if there is an edge in the surface part between the
#     # profile and end of spine. Should just remove all edges up to a point.
#     # Then be sure to use those when creating spines for run-outs.
#
#     # Run-out 1
#     spine2 = ShapeTools.to_wire(ordered_edges[0])
#     shape2 = ShapeTools.make_pipe_shell(spine2, [spine0.v1, profile1],
#                                         surface_part, True)
#
#     # Run-out 2
#     spine3 = ShapeTools.to_wire(ordered_edges[-1])
#     shape3 = ShapeTools.make_pipe_shell(spine3, [profile2, spine0.v2],
#                                         surface_part, True)
#     if not shape2 or not shape3:
#         return None
#
#     # Fuse together shapes.
#     shape = ShapeTools.bop_algo([shape1], [shape2, shape3])
#     if not shape:
#         return None
#
#     # Create stringer.
#     if etype in ['stringer']:
#         part = Stringer(label, shape, cref, sref, assy)
#         return part
#
#     # Create stiffener
#     stiffener = Stiffener2D(label, shape, cref, sref)
#
#     # TODO Fuse the surface part and the stiffener?
#     # surface_part.fuse(stiffener)
#     # bop = ShapeTools.bfuse(surface_part, stiffener, 'builder')
#     # if not bop.IsDone():
#     #     return None
#     # surface_part.set_shape(bop.Shape())
#     #
#     # # Reshape stiffener part.
#     # stiffener.reshape(bop)
#
#     # Add as subpart.
#     # surface_part.add_subpart(stiffener.label, stiffener)
#     # TODO How to handle stiffener 2-d as subpart?
#
#     return stiffener
#
#
# def create_stiffener2d_by_section(etype, label, surface_part, surface_shape, h,
#                                   runout_angle, cut_part=False, assy=None):
#     """
#     Create a 2-D stiffener on a surface part by intersection.
#     """
#     # if not isinstance(surface_part, SurfacePart):
#     #     return None
#
#     if h <= 0. or runout_angle < 0:
#         return None
#
#     sref = None
#     if CheckGeom.is_surface_like(surface_shape):
#         sref = surface_shape
#
#     surface_shape = ShapeTools.to_shape(surface_shape)
#     if not surface_shape:
#         return None
#
#     # Build spine
#     # edges = ShapeTools.bsection(surface_part, surface_shape, 'edge')
#     # if not edges:
#     #     return None
#     bop = ShapeTools.bcut(surface_part, surface_shape, 'builder')
#     if bop.ErrorStatus() != 0:
#         print('here 1')
#         return None
#
#     list_of_edges = bop.SectionEdges()
#     if list_of_edges.IsEmpty():
#         print('here 2')
#         return None
#
#     edges = []
#     while not list_of_edges.IsEmpty():
#         edges.append(list_of_edges.First())
#         list_of_edges.RemoveFirst()
#
#     wires = ShapeTools.connect_edges(edges)
#     if not wires:
#         return None
#     if len(wires) == 1:
#         spine0 = wires[0]
#     else:
#         spine0 = ShapeTools.longest_wire(wires)
#     if not spine0:
#         return None
#
#     stiffener = create_stiffener2d_by_wire(etype, label, surface_part, spine0,
#                                            h, runout_angle, sref, assy)
#     if not stiffener:
#         return None
#
#     # TODO Reshape surface part after stringer cut?
#     if cut_part:
#         surface_part.reshape(bop)
#
#     return stiffener
#
#
# def create_stiffener2d_by_sections(etype, label, surface_part,
#                                    surface_shape, h, runout_angle, assy=None):
#     """
#     Create a 2-D stiffener on a surface part by all intersections.
#     """
#     # if not isinstance(surface_part, SurfacePart):
#     #     return []
#
#     if h <= 0. or runout_angle < 0:
#         return None
#
#     sref = None
#     if CheckGeom.is_surface_like(surface_shape):
#         sref = surface_shape
#
#     surface_shape = ShapeTools.to_shape(surface_shape)
#     if not surface_shape:
#         return []
#
#     # Build spine
#     edges = ShapeTools.bsection(surface_part, surface_shape, 'edge')
#     if not edges:
#         return []
#     wires = ShapeTools.connect_edges(edges)
#     if not wires:
#         return []
#
#     # Create stiffener for all wires.
#     stiffeners = []
#     for wire in wires:
#         stiffener = create_stiffener2d_by_wire(etype, label, surface_part,
#                                                wire, h, runout_angle, sref,
#                                                assy)
#         if not stiffener:
#             continue
#         stiffeners.append(stiffener)
#
#     return stiffeners
