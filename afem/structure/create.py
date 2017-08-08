from OCC.TopoDS import TopoDS_Edge, TopoDS_Face, TopoDS_Shape, TopoDS_Shell, \
    TopoDS_Wire

from afem.geometry.check import CheckGeom
from afem.geometry.create import *
from afem.geometry.entities import *
from afem.oml.entities import Wing
from afem.structure.entities import *
from afem.topology.bop import *
from afem.topology.check import CheckShape
from afem.topology.create import *
from afem.topology.explore import ExploreShape
from afem.topology.props import *

__all__ = ["CurvePartByShape", "BeamByShape", "BeamByCurve", "BeamByPoints",
           "SurfacePartByShape", "SparByParameters", "SparByPoints",
           "SparBySurface", "SparBetweenShapes"]


def _validate_type(input_, type_, required=True):
    """
    Validate the input is of the required type.
    """
    if not required and not input_:
        return None
    if not isinstance(input_, type_):
        msg = 'Invalid input type.'
        raise TypeError(msg)


class CurvePartByShape(object):
    """
    Create a curve part using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> from afem.structure import CurvePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> part = CurvePartByShape('part', e).curve_part
    """

    def __init__(self, label, shape, cref=None):
        _validate_type(label, str)
        _validate_type(shape, (TopoDS_Edge, TopoDS_Wire))
        _validate_type(cref, Curve, False)

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._curve_part = CurvePart(label, shape, cref)

    @property
    def curve_part(self):
        """
        :return: The curve part.
        :rtype: afem.structure.entities.CurvePart
        """
        return self._curve_part


class BeamByShape(object):
    """
    Create a beam using a shape.

    :param str label: Part label.
    :param shape: The shape.
    :type shape: OCC.TopoDS.TopoDS_Edge or OCC.TopoDS.TopoDS_Wire
    :param afem.geometry.entities.Curve cref: The reference curve. If not
        provided then a curve will be extracted from the shape.

    Usage:

    >>> from afem.topology import EdgeByPoints
    >>> from afem.structure import BeamByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> beam = BeamByShape('part', e).beam
    """

    def __init__(self, label, shape, cref=None):
        _validate_type(label, str)
        _validate_type(shape, (TopoDS_Edge, TopoDS_Wire))
        _validate_type(cref, Curve, False)

        if cref is None:
            cref = ExploreShape.curve_of_shape(shape)

        self._beam = Beam(label, shape, cref)

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

    Usage:

    >>> from afem.geometry import NurbsCurveByPoints
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

    >>> from afem.topology import EdgeByPoints, FaceByDrag
    >>> from afem.structure import SurfacePartByShape
    >>> e = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> f = FaceByDrag(e, (0., 10., 0.)).face
    >>> part = SurfacePartByShape('part', f).surface_part
    >>> part.area
    99.99999999999997
    """

    def __init__(self, label, shape, cref=None, sref=None):
        _validate_type(label, str)
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
        :rtype: afem.structure.entities.SurfacePart
        """
        return self._surface_part


class SparByParameters(object):
    """
    Create a spar between wing parameters.

    :param str label: Part label.
    :param float u1: Starting point u-parameter.
    :param float v1: Starting point v-parameter.
    :param float u2: Ending point u-parameter.
    :param float v2: Ending point v-parameter.
    :param afem.oml.entities.Wing wing: The wing.
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
        _validate_type(label, str)
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
    :param afem.oml.entities.Wing wing: The wing.
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
    :param afem.geometry.entities.Surface sref: The reference surface.
    :param afem.oml.entities.Wing wing: The wing.

    .. note::

        The reference curve of the part will be untrimmed in this method. It
        will be determined by the intersection between *sref* and the wing
        reference shape.
    """

    def __init__(self, label, sref, wing):
        _validate_type(label, str)
        _validate_type(sref, Surface)
        _validate_type(wing, Wing)

        basis_shape = FaceBySurface(sref).face

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
        self._spar = Spar(label, shell, cref, sref)

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
    :param OCC.TopoDS.TopoDS_Shape shape1: Starting shape.
    :param OCC.TopoDS.TopoDS_Shape shape2: Ending shape.
    :param afem.oml.entities.Wing wing: The wing.
    :param basis_shape: The basis shape.
    :type basis_shape: OCC.TopoDS.TopoDS_Face or OCC.TopoDS.TopoDS_Shell
    """

    def __init__(self, label, shape1, shape2, wing, basis_shape):
        _validate_type(label, str)
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


class SparsBetweenPlanes(object):
    pass


class SparsAtShapes(object):
    pass


class SparsAlongCurve(object):
    pass


class RibByParameters(object):
    pass


class RibByPoints(object):
    pass


class RibByShape(object):
    pass


class RibBetweenGeom(object):
    pass


class RibsBetweenPlanes(object):
    pass


class RibsAtShapes(object):
    pass


class RibsAlongCurve(object):
    pass


class BulkheadByShape(object):
    pass


class FloorByShape(object):
    pass


class FrameByShape(object):
    pass


class FramesBetweenPlanes(object):
    pass


class FramesAtShapes(object):
    pass


class SkinBySolid(object):
    pass


class SkinByBody(object):
    pass


class StringerBySpine(object):
    pass


class StringerBySection(object):
    pass


class StringersBySections(object):
    pass


if __name__ == "__main__":
    import doctest

    doctest.testmod()
