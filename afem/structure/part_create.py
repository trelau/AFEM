from OCC.TopoDS import TopoDS_Edge, TopoDS_Face, TopoDS_Shape, TopoDS_Shell, \
    TopoDS_Solid, TopoDS_Wire

from afem.geometry.geom_check import CheckGeom
from afem.geometry.geom_create import *
from afem.geometry.geom_entities import *
from afem.oml.oml_entities import Body, Fuselage, Wing
from afem.structure.part_entities import *
from afem.topology.topo_bop import *
from afem.topology.topo_check import CheckShape
from afem.topology.topo_create import *
from afem.topology.topo_explore import ExploreFreeEdges, ExploreShape
from afem.topology.topo_props import *

__all__ = ["CurvePartByShape", "BeamByShape", "BeamByCurve", "BeamByPoints",
           "SurfacePartByShape", "SparByParameters", "SparByPoints",
           "SparBySurface", "SparBetweenShapes", "RibByParameters",
           "RibByPoints", "RibBySurface", "RibBetweenShapes",
           "RibsBetweenPlanesByNumber", "RibsBetweenPlanesByDistance",
           "RibsAlongCurveByNumber", "RibsAlongCurveByDistance",
           "RibsAtShapes", "BulkheadBySurface", "FloorBySurface",
           "FrameByPlane", "FramesBetweenPlanesByNumber",
           "FramesBetweenPlanesByDistance", "FramesByPlanes", "SkinBySolid",
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
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.

    Usage:

    >>> from afem.topo_patch import EdgeByPoints
    >>> from afem.structure import CurvePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> part = CurvePartByShape('part', e).curve_part
    """

    def __init__(self, label, shape, cref=None):
        _validate_type(shape, (TopoDS_Edge, TopoDS_Wire))
        _validate_type(cref, Curve, False)

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._curve_part = CurvePart(label, shape, cref)

    @property
    def curve_part(self):
        """
        :return: The curve part.
        :rtype: afem.structure.part_entities.CurvePart
        """
        return self._curve_part


# BEAM -----------------------------------------------------------------------

class BeamByShape(object):
    """
    Create a beam using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.

    Usage:

    >>> from afem.topo_patch import EdgeByPoints
    >>> from afem.structure import BeamByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> beam = BeamByShape('part', e).beam
    """

    def __init__(self, label, shape, cref=None):
        _validate_type(shape, (TopoDS_Edge, TopoDS_Wire))
        _validate_type(cref, Curve, False)

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._beam = Beam(label, shape, cref)

    @property
    def beam(self):
        """
        :return: The beam.
        :rtype: afem.structure.part_entities.Beam
        """
        return self._beam


class BeamByCurve(BeamByShape):
    """
    Create a beam using a curve.

    :param str label: Part label.
    :param afem.geometry.entities.Curve crv: The curve.

    Usage:

    >>> from afem.geom_patch import NurbsCurveByPoints
    >>> from afem.structure import BeamByCurve
    >>> c = NurbsCurveByPoints([(0., 0., 0.), (10., 0., 0.)]).curve
    >>> part = BeamByCurve('part', c).beam
    """

    def __init__(self, label, crv):
        _validate_type(crv, Curve)

        e = EdgeByCurve(crv).edge
        super(BeamByCurve, self).__init__(label, e, crv)


class BeamByPoints(BeamByCurve):
    """
    Create a beam between two points.

    :param str label: Part label.
    :param point_like p1: First point.
    :param point_like p2: Second point.

    Usage:

    >>> from afem.structure import BeamByPoints
    >>> beam = BeamByPoints('part', (0., 0., 0.), (10., 0., 0.)).beam
    >>> beam.length
    10.0
    """

    def __init__(self, label, p1, p2):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        _validate_type(p1, Point)
        _validate_type(p2, Point)

        c = NurbsCurveByPoints([p1, p2]).curve
        super(BeamByPoints, self).__init__(label, c)


# SURFACE PART ----------------------------------------------------------------

class SurfacePartByShape(object):
    """
    Create a surface part using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    :param afem.geometry.entities.Curve cref: The reference curve.
    :param afem.geometry.entities.Surface sref: The reference surface. If
        not provided then the basis surface of the largest face of the shape
        will be used.

    Usage:

    >>> from afem.topo_patch import EdgeByPoints, FaceByDrag
    >>> from afem.structure import SurfacePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> f = FaceByDrag(e, (0., 10., 0.)).face
    >>> part = SurfacePartByShape('part', f).surface_part
    >>> part.area
    99.99999999999997
    """

    def __init__(self, label, shape, cref=None, sref=None):
        _validate_type(shape, (TopoDS_Face, TopoDS_Shell))
        _validate_type(cref, Curve, False)
        _validate_type(sref, Surface, False)

        if sref is None:
            sref = ExploreShape.surface_of_shape(shape)

        self._surface_part = SurfacePart(label, shape, cref, sref)

    @property
    def surface_part(self):
        """
        :return: The surface part.
        :rtype: afem.structure.part_entities.SurfacePart
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
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.

    .. note::

        * If a basis shape is provided, then the starting and ending parameters
          should be near the intersection between this shape and the wing
          reference shape.
    """

    def __init__(self, label, u1, v1, u2, v2, wing, basis_shape=None):
        _validate_type(wing, Wing)
        _validate_type(basis_shape, (TopoDS_Face, TopoDS_Shell), False)

        # Determine reference surface and basis shape
        if basis_shape is None:
            pln = wing.extract_plane(u1, v1, u2, v2)
            sref = pln
            basis_shape = FaceBySurface(pln).face
        else:
            sref = ExploreShape.surface_of_shape(basis_shape)

        # Extract wing cref
        cref = wing.extract_curve(u1, v1, u2, v2, basis_shape)

        # Build part shape
        common = CommonShapes(basis_shape, wing)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape
        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._spar = Spar(label, shell, cref, sref)

    @property
    def spar(self):
        """
        :return: The spar.
        :rtype: afem.structure.part_entities.Spar
        """
        return self._spar


class SparByPoints(SparByParameters):
    """
    Create a spar between two points. This method inverts the starting and
    ending points and then uses :class:`.SparByParameters`.

    :param str label: Part label.
    :param point_like p1: Starting point.
    :param point_like p2: End point.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell

    .. note::

        * The starting and ending points should be on or near the wing
          reference surface since they are projected to find parameters.

        * If a basis shape is provided, then the starting and ending points
          should be near the intersection between this shape and the wing
          reference shape.
    """

    def __init__(self, label, p1, p2, wing, basis_shape=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        _validate_type(p1, Point)
        _validate_type(p2, Point)

        # Invert points
        u1, v1 = wing.invert(p1)
        u2, v2 = wing.invert(p2)

        # Use SparByParameters
        super(SparByPoints, self).__init__(label, u1, v1, u2, v2, wing,
                                           basis_shape)


class SparBySurface(object):
    """
    Create a spar using a surface.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Surface srf: The surface.
    :param afem.oml.oml_entities.Wing wing: The wing.

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *srf* and the wing
        reference shape.
    """

    def __init__(self, label, srf, wing):
        _validate_type(srf, Surface)
        _validate_type(wing, Wing)

        basis_shape = FaceBySurface(srf).face

        # Build reference curve
        section = IntersectShapes(basis_shape, wing.sref_shape)
        edges = ExploreShape.get_edges(section.shape)
        wires = WiresByConnectedEdges(edges).wires
        w = LengthOfShapes(wires).longest_shape
        w = CheckShape.to_wire(w)
        cref = ExploreShape.curve_of_shape(w)

        # Orient cref so that p1 is nearest (u1, v1) on the wing
        umin, vmin = wing.sref.u1, wing.sref.v1
        p0 = wing.eval(umin, vmin)
        p1 = cref.eval(cref.u1)
        p2 = cref.eval(cref.u2)
        d1 = p0.distance(p1)
        d2 = p0.distance(p2)
        if d2 < d1:
            cref.reverse()

        # Build part shape
        common = CommonShapes(basis_shape, wing)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape
        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._spar = Spar(label, shell, cref, srf)

    @property
    def spar(self):
        """
        :return: The spar.
        :rtype: afem.structure.part_entities.Spar
        """
        return self._spar


class SparBetweenShapes(SparByPoints):
    """
    Create a spar between shapes. This method uses the shapes to define
    starting and ending points by intersecting the reference curve and then
    uses :class:`.SparByPoints`.

    :param str label: Part label.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param basis_shape: The basis shape.
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    """

    def __init__(self, label, shape1, shape2, wing, basis_shape):
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(wing, Wing)
        _validate_type(basis_shape, (TopoDS_Face, TopoDS_Shell), False)

        wing_basis_edges = IntersectShapes(basis_shape, wing.sref_shape).shape
        p1_shape = IntersectShapes(shape1, wing_basis_edges).shape
        p2_shape = IntersectShapes(shape2, wing_basis_edges).shape
        v1 = ExploreShape.get_vertices(p1_shape)[0]
        v2 = ExploreShape.get_vertices(p2_shape)[0]
        p1 = ExploreShape.pnt_of_vertex(v1)
        p2 = ExploreShape.pnt_of_vertex(v2)

        super(SparBetweenShapes, self).__init__(label, p1, p2, wing,
                                                basis_shape)


class SparsBetweenPlanesByNumber(object):
    # TODO SparsBetweenPlanesByNumber
    pass


class SparsBetweenPlanesByDistance(object):
    # TODO SparsBetweenPlanesByDistance
    pass


class SparsAlongCurveByNumber(object):
    # TODO SparsAlongCurveByNumber
    pass


class SparsAlongCurveByDistance(object):
    # TODO SparsAlongCurveByDistance
    pass


class SparsAtShapes(object):
    # TODO SparsAtShapes
    pass


# RIB -------------------------------------------------------------------------

class RibByParameters(object):
    """
    Create a rib between wing parameters.

    :param str label: Part label.
    :param float u1: Starting point u-parameter.
    :param float v1: Starting point v-parameter.
    :param float u2: Ending point u-parameter.
    :param float v2: Ending point v-parameter.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell

    .. note::

        * If a basis shape is provided, then the starting and ending parameters
          should be near the intersection between this shape and the wing
          reference shape.
    """

    def __init__(self, label, u1, v1, u2, v2, wing, basis_shape=None):
        _validate_type(wing, Wing)
        _validate_type(basis_shape, (TopoDS_Face, TopoDS_Shell), False)

        # Determine reference surface and basis shape
        if basis_shape is None:
            pln = wing.extract_plane(u1, v1, u2, v2)
            sref = pln
            basis_shape = FaceBySurface(pln).face
        else:
            sref = ExploreShape.surface_of_shape(basis_shape)

        # Extract wing cref
        cref = wing.extract_curve(u1, v1, u2, v2, basis_shape)

        # Build part shape
        common = CommonShapes(basis_shape, wing)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape
        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._rib = Rib(label, shell, cref, sref)

    @property
    def rib(self):
        """
        :return: The rib.
        :rtype: afem.structure.part_entities.Rib
        """
        return self._rib


class RibByPoints(RibByParameters):
    """
    Create a rib between two points. This method inverts the starting and
    ending points and then uses :class:`.RibByParameters`.

    :param str label: Part label.
    :param point_like p1: Starting point.
    :param point_like p2: End point.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param basis_shape: The basis shape to define the shape of the part. If
        not provided, then a plane will be defined between (u1, v1),
        (u2, v2), and a point translated from the reference surface normal at
        (u1, v1).
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell

    .. note::

        * The starting and ending points should be on or near the wing
          reference surface since they are projected to find parameters.

        * If a basis shape is provided, then the starting and ending points
          should be near the intersection between this shape and the wing
          reference shape.
    """

    def __init__(self, label, p1, p2, wing, basis_shape=None):
        p1 = CheckGeom.to_point(p1)
        p2 = CheckGeom.to_point(p2)

        _validate_type(p1, Point)
        _validate_type(p2, Point)

        # Invert points
        u1, v1 = wing.invert(p1)
        u2, v2 = wing.invert(p2)

        # Use SparByParameters
        super(RibByPoints, self).__init__(label, u1, v1, u2, v2, wing,
                                          basis_shape)


class RibBySurface(object):
    """
    Create a rib using a surface.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Surface srf: The surface.
    :param afem.oml.oml_entities.Wing wing: The wing.

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *srf* and the wing
        reference shape.
    """

    def __init__(self, label, srf, wing):
        _validate_type(srf, Surface)
        _validate_type(wing, Wing)

        basis_shape = FaceBySurface(srf).face

        # Build reference curve
        section = IntersectShapes(basis_shape, wing.sref_shape)
        edges = ExploreShape.get_edges(section.shape)
        wires = WiresByConnectedEdges(edges).wires
        w = LengthOfShapes(wires).longest_shape
        w = CheckShape.to_wire(w)
        cref = ExploreShape.curve_of_shape(w)

        # Orient cref so that p1 is nearest (u1, v1) on the wing
        umin, vmin = wing.sref.u1, wing.sref.v1
        p0 = wing.eval(umin, vmin)
        p1 = cref.eval(cref.u1)
        p2 = cref.eval(cref.u2)
        d1 = p0.distance(p1)
        d2 = p0.distance(p2)
        if d2 < d1:
            cref.reverse()

        # Build part shape
        common = CommonShapes(basis_shape, wing)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape
        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._rib = Rib(label, shell, cref, srf)

    @property
    def rib(self):
        """
        :return: The rib.
        :rtype: afem.structure.part_entities.Rib
        """
        return self._rib


class RibBetweenShapes(RibByPoints):
    """
    Create a rib between shapes. This method uses the shapes to define
    starting and ending points by intersecting the reference curve and then
    uses :class:`.RibByPoints`.

    :param str label: Part label.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param basis_shape: The basis shape.
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    """

    def __init__(self, label, shape1, shape2, wing, basis_shape):
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(wing, Wing)
        _validate_type(basis_shape, (TopoDS_Face, TopoDS_Shell), False)

        wing_basis_edges = IntersectShapes(basis_shape, wing.sref_shape).shape
        p1_shape = IntersectShapes(shape1, wing_basis_edges).shape
        p2_shape = IntersectShapes(shape2, wing_basis_edges).shape
        v1 = ExploreShape.get_vertices(p1_shape)[0]
        v2 = ExploreShape.get_vertices(p2_shape)[0]
        p1 = ExploreShape.pnt_of_vertex(v1)
        p2 = ExploreShape.pnt_of_vertex(v2)

        super(RibBetweenShapes, self).__init__(label, p1, p2, wing,
                                               basis_shape)


class RibsBetweenPlanesByNumber(object):
    """
    Create a specified number of planar ribs between two planes. This method
    uses :class:`.PlanesBetweenPlanesByNumber` and :class:`.RibBetweenShapes`.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Plane pln1: The first plane.
    :param afem.geometry.geom_entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param float d1: An offset distance for the first plane. This is typically
        a positive number indicating a distance from *u1* towards *u2*.
    :param float d2: An offset distance for the last plane. This is typically
        a negative number indicating a distance from *u2* towards *u1*.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    """

    def __init__(self, label, pln1, pln2, n, shape1, shape2, wing, d1=None,
                 d2=None, first_index=1, delimiter=' '):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(wing, Wing)
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
            rib = RibBetweenShapes(label_indx, shape1, shape2, wing,
                                   basis_shape).rib
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
        :rtype: list[afem.structure.part_entities.Rib]
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
    :param afem.geometry.geom_entities.Plane pln1: The first plane.
    :param afem.geometry.geom_entities.Plane pln2: The second plane.
    :param float maxd: The maximum allowed spacing. The actual spacing will
        be adjusted to not to exceed this value.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.oml_entities.Wing wing: The wing.
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
    """

    def __init__(self, label, pln1, pln2, maxd, shape1, shape2, wing, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' '):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(wing, Wing)
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
            rib = RibBetweenShapes(label_indx, shape1, shape2, wing,
                                   basis_shape).rib
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
        :rtype: list[afem.structure.part_entities.Rib]
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
    :param afem.geometry.geom_entities.Curve crv: The curve.
    :param int n: The number of parts.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param afem.geometry.geom_entities.Plane ref_pln: The normal of this plane
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
    """

    def __init__(self, label, crv, n, shape1, shape2, wing, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, first_index=1,
                 delimiter=' ', tol=1.0e-7):
        _validate_type(crv, Curve)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(wing, Wing)
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
            rib = RibBetweenShapes(label_indx, shape1, shape2, wing,
                                   basis_shape).rib
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
        :rtype: list[afem.structure.part_entities.Rib]
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
    :param afem.geometry.geom_entities.Curve crv: The curve.
    :param float maxd: The maximum allowed spacing between planes. The
        actual spacing will be adjusted to not to exceed this value.
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.oml_entities.Wing wing: The wing.
    :param afem.geometry.geom_entities.Plane ref_pln: The normal of this plane
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
    """

    def __init__(self, label, crv, maxd, shape1, shape2, wing, ref_pln=None,
                 u1=None, u2=None, d1=None, d2=None, nmin=0, first_index=1,
                 delimiter=' ', tol=1.0e-7):
        _validate_type(crv, Curve)
        _validate_type(shape1, TopoDS_Shape)
        _validate_type(shape2, TopoDS_Shape)
        _validate_type(wing, Wing)
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
            rib = RibBetweenShapes(label_indx, shape1, shape2, wing,
                                   basis_shape).rib
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
        :rtype: list[afem.structure.part_entities.Rib]
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


class RibsAtShapes(object):
    # TODO RibsAtShapes
    pass


# BULKHEAD -------------------------------------------------------------------

class BulkheadBySurface(object):
    """
    Create a bulkhead using a surface.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Surface srf: The surface.
    :param afem.oml.oml_entities.Fuselage fuselage: The fuselage.

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, srf, fuselage):
        _validate_type(srf, Surface)
        _validate_type(fuselage, Fuselage)

        basis_shape = FaceBySurface(srf).face

        # Build part shape
        common = CommonShapes(basis_shape, fuselage)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape
        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._bh = Bulkhead(label, shell, sref=srf)

    @property
    def bulkhead(self):
        """
        :return: The bulkhead.
        :rtype: afem.structure.part_entities.Bulkhead
        """
        return self._bh


# FLOOR -----------------------------------------------------------------------

class FloorBySurface(object):
    """
    Create a floor using a surface.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Surface srf: The surface.
    :param afem.oml.oml_entities.Fuselage fuselage: The fuselage.

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, srf, fuselage):
        _validate_type(srf, Surface)
        _validate_type(fuselage, Fuselage)

        basis_shape = FaceBySurface(srf).face

        # Build part shape
        common = CommonShapes(basis_shape, fuselage)
        if not common.is_done:
            msg = 'Boolean operation failed.'
            raise RuntimeError(msg)
        shape = common.shape
        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._floor = Floor(label, shell, sref=srf)

    @property
    def floor(self):
        """
        :return: The floor.
        :rtype: afem.structure.part_entities.Floor
        """
        return self._floor


# FRAME -----------------------------------------------------------------------

class FrameByPlane(object):
    """
    Create a frame using a plane. A plane is required since the shape of
    frame is formed using an offset algorithm that requires planar shapes.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Plane pln: The plane.
    :param afem.oml.oml_entities.Fuselage fuselage: The fuselage.
    :param float height: The height. The absolute value is used.

    :raise TypeError: If an input type is invalid.
    :raise RuntimeError: If Boolean operation failed.
    """

    def __init__(self, label, pln, fuselage, height):
        _validate_type(pln, Plane)
        _validate_type(fuselage, Fuselage)

        basis_shape = FaceBySurface(pln).face

        # Find initial shape
        common = CommonShapes(basis_shape, fuselage)
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

        faces = ExploreShape.get_faces(shape)
        shell = ShellByFaces(faces).shell

        # Create the part
        self._frame = Frame(label, shell, sref=pln)

    @property
    def frame(self):
        """
        :return: The frame.
        :rtype: afem.structure.part_entities.Frame
        """
        return self._frame


class FramesBetweenPlanesByNumber(object):
    """
    Create a specified number of frames between two planes. This method uses
    :class:`.PlanesBetweenPlanesByNumber` and :class:`.FrameByPlane`.

    :param str label: Part label.
    :param afem.geometry.geom_entities.Plane pln1: The first plane.
    :param afem.geometry.geom_entities.Plane pln2: The second plane.
    :param int n: The number of parts.
    :param afem.oml.oml_entities.Fuselage fuselage: The fuselage.
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
    """

    def __init__(self, label, pln1, pln2, n, fuselage, height, d1=None,
                 d2=None, first_index=1, delimiter=' '):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(fuselage, Fuselage)

        n = int(n)
        first_index = int(first_index)

        builder = PlanesBetweenPlanesByNumber(pln1, pln2, n, d1, d2)

        self._frames = []
        self._nframes = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, fuselage, height).frame
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
        :rtype: list[afem.structure.part_entities.Frame]
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
    :param afem.geometry.geom_entities.Plane pln1: The first plane.
    :param afem.geometry.geom_entities.Plane pln2: The second plane.
    :param float maxd: The maximum allowed spacing. The actual spacing will
        be adjusted to not to exceed this value.
    :param afem.oml.oml_entities.Fuselage fuselage: The fuselage.
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
    """

    def __init__(self, label, pln1, pln2, maxd, fuselage, height, d1=None,
                 d2=None, nmin=0, first_index=1, delimiter=' '):
        _validate_type(pln1, Plane)
        _validate_type(pln2, Plane)
        _validate_type(fuselage, Fuselage)

        first_index = int(first_index)

        builder = PlanesBetweenPlanesByDistance(pln1, pln2, maxd, d1, d2, nmin)

        self._frames = []
        self._nframes = builder.nplanes
        self._ds = builder.spacing
        for pln in builder.planes:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, fuselage, height).frame
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
        :rtype: list[afem.structure.part_entities.Frame]
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


class FramesByPlanes(object):
    """
    Create frames using a list of planes. This method uses
    :class:`.FrameByPlane`.

    :param str label: Part label.
    :param list[afem.geometry.geom_entities.Plane] plns: The planes.
    :param afem.oml.oml_entities.Fuselage fuselage: The fuselage.
    :param float height: The height.
    :param int first_index: The first index appended to the part label as
        parts are created successively.
    :param str delimiter: The delimiter to use when joining the part label
        with the index. The final part label will be
        'label' + 'delimiter' + 'index'.
    """

    def __init__(self, label, plns, fuselage, height, first_index=1,
                 delimiter=' '):
        for pln in plns:
            _validate_type(pln, Plane)
        _validate_type(fuselage, Fuselage)

        first_index = int(first_index)

        self._frames = []
        self._nframes = len(plns)
        self._ds = None
        for pln in plns:
            label_indx = delimiter.join([label, str(first_index)])
            frame = FrameByPlane(label_indx, pln, fuselage, height).frame
            first_index += 1
            self._frames.append(frame)
        self._next_index = first_index

        if len(plns) > 1:
            p1 = plns[0].eval(0., 0.)
            p2 = plns[1].eval(0., 0.)
            self._ds = p1.distance(p2)

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
        :rtype: list[afem.structure.part_entities.Frame]
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
    """

    def __init__(self, label, solid, copy=True):
        _validate_type(solid, TopoDS_Solid)

        shell = ExploreShape.outer_shell(solid)
        if copy:
            shell = ExploreShape.copy_shape(shell, False)

        self._skin = Skin(label, shell)

    @property
    def skin(self):
        """
        :return: The skin.
        :rtype: afem.structure.part_entities.Skin
        """
        return self._skin


class SkinByBody(SkinBySolid):
    """
    Create a skin part from the outer shell of a body.

    :param str label: Part label.
    :param afem.oml.oml_entities.Body body: The body.
    :param bool copy: Option to copy the outer shell.
    """

    def __init__(self, label, body, copy=True):
        _validate_type(body, Body)
        super(SkinByBody, self).__init__(label, body, copy)


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
