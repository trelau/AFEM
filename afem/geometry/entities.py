from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.BRepGProp import brepgprop_SurfaceProperties
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.GProp import GProp_GProps
from OCC.Geom import (Handle_Geom_BSplineCurve, Handle_Geom_BSplineSurface,
                      Handle_Geom_Circle, Handle_Geom_Ellipse,
                      Handle_Geom_Line,
                      Handle_Geom_Plane, Handle_Geom_TrimmedCurve)
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface
from OCC.TColStd import (TColStd_Array1OfInteger, TColStd_Array1OfReal,
                         TColStd_Array2OfReal)
from OCC.TColgp import (TColgp_Array1OfPnt, TColgp_Array2OfPnt)
from OCC.gp import gp_Ax1, gp_Ax3, gp_Dir, gp_Pnt, gp_Pnt2d, gp_Vec, gp_XYZ
from numpy import add, array, float64, ndarray, subtract

from afem.geometry.utils import (global_to_local_param,
                                 homogenize_array1d,
                                 homogenize_array2d,
                                 local_to_global_param,
                                 reparameterize_knots)
from afem.graphics.viewer import ViewableItem
from afem.misc.check import is_array_like
from afem.occ.utils import (to_np_from_tcolgp_array1_pnt,
                            to_np_from_tcolgp_array2_pnt,
                            to_np_from_tcolstd_array1_integer,
                            to_np_from_tcolstd_array1_real,
                            to_np_from_tcolstd_array2_real,
                            to_tcolgp_array1_pnt, to_tcolstd_array1_real)

__all__ = ["Geometry2D", "Point2D", "NurbsCurve2D", "Geometry", "Point",
           "Direction", "Vector", "Axis1", "Axis3", "Curve", "Line", "Circle",
           "NurbsCurve", "TrimmedCurve", "Surface", "Plane", "NurbsSurface"]


# 2-D -------------------------------------------------------------------------

class Geometry2D(ViewableItem):
    """
    Base class for 2-D geometry.

    :param OCC.Geom2d.Handle_Geom2d_Geometry the_handle: The geometry handle.
    """

    def __init__(self, the_handle=None):
        super(Geometry2D, self).__init__()
        self._handle = the_handle
        self._object = None
        if the_handle is not None:
            self._object = the_handle.GetObject()


class Point2D(gp_Pnt2d, Geometry2D):
    """
    A 2-D Cartesian point. Supports NumPy array methods.

    For more information see gp_Pnt2d_.

    .. _gp_Pnt2d: https://www.opencascade.com/doc/occt-7.1.0/refman/htmlclassgp___pnt2d.html

    Usage:

    >>> from afem.geometry import Point2D
    >>> Point2D()
    Point2D(0.0, 0.0)
    >>> Point2D(1., 2.)
    Point2D(1.0, 2.0)
    >>> from numpy import array
    >>> array(Point2D(1., 2.))
    array([ 1.,  2.])
    >>> p1 = Point2D(1., 2.)
    >>> p2 = Point2D(4., 5.)
    >>> p1[0]
    1.0
    >>> p1[1]
    2.0
    >>> p1 + p2
    array([ 5.,  7.])
    >>> p2 - p1
    array([ 3.,  3.])
    >>> array([p1, p2])
    array([[ 1.,  2.],
           [ 4.,  5.]])
    """

    def __init__(self, *args):
        super(Point2D, self).__init__(*args)
        Geometry2D.__init__(self)

    def __str__(self):
        return 'Point2D({0}, {1})'.format(*self.xy)

    def __repr__(self):
        return 'Point2D({0}, {1})'.format(*self.xy)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.xy, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.xy:
            yield elm

    def __len__(self):
        return 2

    def __getitem__(self, item):
        return self.xy[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def xy(self):
        """
        :return: The point xy-location.
        :rtype: numpy.ndarray
        """
        return array([self.X(), self.Y()], dtype=float64)

    @property
    def x(self):
        return self.X()

    @x.setter
    def x(self, x):
        """
        The point x-location.

        :getter: Returns the x-location.
        :setter: Sets the x-location.
        :type: float
        """
        self.SetX(x)

    @property
    def y(self):
        """
        The point y-location.

        :getter: Returns the y-location.
        :setter: Sets the y-location.
        :type: float
        """
        return self.Y()

    @y.setter
    def y(self, y):
        self.SetY(y)

    def set_xy(self, xy):
        """
        Set point coordinates.

        :param array_like xy: Point coordinates.

        :return: *True* if set, *False* if not.
        :rtype: bool
        """
        if isinstance(xy, gp_Pnt):
            self.SetXY(xy.XYZ())
            return True
        if isinstance(xy, gp_XYZ):
            self.SetXY(xy)
            return True
        if is_array_like(xy) and len(xy) == 2:
            self.x, self.y = xy
            return True
        return False

    def distance(self, other):
        """
        Compute the distance between two points.

        :param point_like other: The other point.

        :return: Distance to the other point.
        :rtype: float
        """
        if isinstance(other, gp_Pnt2d):
            return self.Distance(other)
        if is_array_like(other) and len(other) == 2:
            other = Point2D(*other)
            return self.Distance(other)
        return None

    def is_equal(self, other, tol=1.0e-7):
        """
        Check for coincident points.

        :param point_like other: The other point.
        :param float tol: Tolerance for coincidence.

        :return: *True* if coincident, *False* if not.
        :rtype: bool
        """
        if isinstance(other, gp_Pnt2d):
            return self.IsEqual(other, tol)
        if is_array_like(other) and len(other) == 2:
            other = Point2D(*other)
            return self.IsEqual(other, tol)
        return False

    def copy(self):
        """
        Return a new copy of the point.

        :return: New point.
        :rtype: afem.geometry.entities.Point2D
        """
        return Point2D(*self.xy)


class NurbsCurve2D(Geometry2D):
    """
    NURBS curve in 2-D space.

    :param OCC.Geom2d.Handle_Geom2d_BSplineCurve the_handle: The curve handle.

    For more information see Geom2d_BSplineCurve_.

    .. _Geom2d_BSplineCurve: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom2d___b_spline_curve.html
    """

    def __init__(self, the_handle):
        super(NurbsCurve2D, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom2d.Handle_Geom2d_BSplineCurve
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom2d.Geom2d_BSplineCurve
        """
        return self._object

    @property
    def adaptor(self):
        """
        :return: A curve adaptor.
        :rtype: OCC.Geom2dAdaptor.Geom2dAdaptor_Curve
        """
        return Geom2dAdaptor_Curve(self.handle)

    @property
    def u1(self):
        """
        :return: The first parameter.
        :rtype: float
        """
        return self.object.FirstParameter()

    @property
    def u2(self):
        """
        :return: The last parameter.
        :rtype: float
        """
        return self.object.LastParameter()

    @property
    def p(self):
        """
        :return: Degree of curve.
        :rtype: int
        """
        return self.object.Degree()

    @property
    def n(self):
        """
        :return: Number of control points - 1.
        :rtype: int
        """
        return self.object.NbPoles() - 1

    @property
    def knots(self):
        """
        :return: Knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.object.NbKnots())
        self.object.Knots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def mult(self):
        """
        :return: Multiplicity of knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.object.NbKnots())
        self.object.Multiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        """
        :return: Knot sequence.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.object.NbPoles() +
                                             self.object.Degree() + 1)
        self.object.KnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def is_closed(self):
        """
        :return: *True* if curve is closed, *False* if not.
        :rtype: bool
        """
        return self.object.IsClosed()

    @property
    def is_periodic(self):
        """
        :return: *True* if curve is periodic, *False* if not.
        :rtype: bool
        """
        return self.object.IsPeriodic()

    def set_domain(self, u1=0., u2=1.):
        """
        Reparameterize the knot vector between *u1* and *u2*.

        :param float u1: First parameter.
        :param float u2: Last parameter.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if u1 > u2:
            return False
        tcol_knots = TColStd_Array1OfReal(1, self.object.NbKnots())
        self.object.Knots(tcol_knots)
        reparameterize_knots(u1, u2, tcol_knots)
        self.object.SetKnots(tcol_knots)
        return True

    def local_to_global_param(self, *args):
        """
        Convert parameter(s) from local domain 0. <= u <= 1. to global domain
        a <= u <= b.

        :param float args: Local parameter(s).

        :return: Global parameter(s).
        :rtype: float or list[float]
        """
        return local_to_global_param(self.u1, self.u2, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0. <= u <= 1.

        :param float args: Global parameter(s).

        :return: Local parameter(s).
        :rtype: float or list[float]
        """
        return global_to_local_param(self.u1, self.u2, *args)

    def eval(self, u):
        """
        Evaluate a point on the curve.

        :param float u: Curve parameter.

        :return: Curve point.
        :rtype: afem.geometry.entities.Point2D
        """
        p = Point2D()
        self.object.D0(u, p)
        return p

    def reverse(self):
        """
        Reverse curve direction.

        :return: None.
        """
        self.object.Reverse()

    def reversed_u(self, u):
        """
        Calculate the parameter on the reversed curve.

        :param float u: Curve parameter.

        :return: Reversed parameter.
        :rtype: float
        """
        return self.object.ReversedParameter(u)

    def segment(self, u1, u2):
        """
        Segment the curve between parameters.

        :param float u1: First parameter.
        :param float u2: Last parameter.

        :return: *True* if segmented, *False* if not.
        :rtype: bool
        """
        if u1 > u2:
            return False
        self.object.Segment(u1, u2)
        return True

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.NurbsCurve2D
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return NurbsCurve2D(h_crv)


# 3-D -------------------------------------------------------------------------

class Geometry(ViewableItem):
    """
    Base class for geometry.

    :param OCC.Geom.Handle_Geom_Geometry the_handle: The geometry handle.
    """

    def __init__(self, the_handle=None):
        super(Geometry, self).__init__()
        self._handle = the_handle
        self._object = None
        if the_handle is not None:
            self._object = the_handle.GetObject()


class Point(gp_Pnt, Geometry):
    """
    A 3-D Cartesian point. Supports NumPy array methods.

    For more information see gp_Pnt_.

    .. _gp_Pnt: https://www.opencascade.com/doc/occt-7.1.0/refman/html/classgp___pnt.html

    Usage:

    >>> from afem.geometry import Point
    >>> Point()
    Point(0.0, 0.0, 0.0)
    >>> Point(1., 2., 3.)
    Point(1.0, 2.0, 3.0)
    >>> from numpy import array
    >>> array(Point(1., 2., 3.))
    array([ 1.,  2.,  3.])
    >>> p1 = Point(1., 2., 3)
    >>> p2 = Point(4., 5., 6.)
    >>> p1[0]
    1.0
    >>> p1[1]
    2.0
    >>> p1[2]
    3.0
    >>> p1 + p2
    array([ 5.,  7.,  9.])
    >>> p2 - p1
    array([ 3.,  3.,  3.])
    >>> array([p1, p2])
    array([[ 1.,  2.,  3.],
           [ 4.,  5.,  6.]])
    """

    def __init__(self, *args):
        super(Point, self).__init__(*args)
        Geometry.__init__(self)

    def __str__(self):
        return 'Point({0}, {1}, {2})'.format(*self.xyz)

    def __repr__(self):
        return 'Point({0}, {1}, {2})'.format(*self.xyz)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.xyz, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.xyz:
            yield elm

    def __len__(self):
        return 3

    def __getitem__(self, item):
        return self.xyz[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def x(self):
        """
        The point x-location.

        :getter: Returns the x-location.
        :setter: Sets the x-location.
        :type: float
        """
        return self.X()

    @x.setter
    def x(self, x):
        self.SetX(x)

    @property
    def y(self):
        """
        The point y-location.

        :getter: Returns the y-location.
        :setter: Sets the y-location.
        :type: float
        """
        return self.Y()

    @y.setter
    def y(self, y):
        self.SetY(y)

    @property
    def z(self):
        """
        The point z-location.

        :getter: Returns the z-location.
        :setter: Sets the z-location.
        :type: float
        """
        return self.Z()

    @z.setter
    def z(self, z):
        self.SetZ(z)

    @property
    def xyz(self):
        """
        :return: The point xyz-location.
        :rtype: numpy.ndarray
        """
        return array([self.X(), self.Y(), self.Z()], dtype=float64)

    def set_xyz(self, xyz):
        """
        Set point coordinates.

        :param point_like xyz: Point coordinates.

        :return: *True* if set, *False* if not.
        :rtype: bool
        """
        if isinstance(xyz, gp_Pnt):
            self.SetXYZ(xyz.XYZ())
            return True
        if isinstance(xyz, gp_XYZ):
            self.SetXYZ(xyz)
            return True
        if is_array_like(xyz) and len(xyz) == 3:
            self.x, self.y, self.z = xyz
            return True
        return False

    def distance(self, other):
        """
        Compute the distance between two points.

        :param point_like other: The other point.

        :return: Distance to the other point.
        :rtype: float
        """
        if isinstance(other, gp_Pnt):
            return self.Distance(other)
        if is_array_like(other) and len(other) == 3:
            other = Point(*other)
            return self.Distance(other)
        return None

    def is_equal(self, other, tol=1.0e-7):
        """
        Check for coincident points.

        :param point_like other: The other point.
        :param float tol: Tolerance for coincidence.

        :return: *True* if coincident, *False* if not.
        :rtype: bool
        """
        if isinstance(other, gp_Pnt):
            return self.IsEqual(other, tol)
        if is_array_like(other) and len(other) == 3:
            other = Point(*other)
            return self.IsEqual(other, tol)
        return False

    def translate(self, v):
        """
        Translate the point along the vector.

        :param vector_like v: The translation vector.

        :return: *True* if translated, *False* if not.
        :rtype: bool

        :raise TypeError: If *v* cannot be converted to a vector.
        """
        if isinstance(v, (tuple, list, ndarray)):
            v = Vector(v[0], v[1], v[2])
        elif isinstance(v, Direction):
            v = Vector(v)
        elif isinstance(v, gp_Vec):
            v = Vector(v.XYZ())
        else:
            msg = 'Invalid vector type.'
            raise TypeError(msg)

        self.Translate(v)
        return True

    def copy(self):
        """
        Return a new copy of the point.

        :return: New point.
        :rtype: afem.geometry.entities.Point
        """
        return Point(*self.xyz)


class Direction(gp_Dir, Geometry):
    """
    Unit vector in 3-D space. Supports NumPy array methods.

    For more information see gp_Dir_.

    .. _gp_Dir: https://www.opencascade.com/doc/occt-7.1.0/refman/html/classgp___dir.html

    Usage:

    >>> from afem.geometry import Direction, Vector
    >>> Direction(10., 0., 0.)
    Direction(1.0, 0.0, 0.0)
    >>> v = Vector(10., 0., 0.)
    >>> Direction(v)
    Direction(1.0, 0.0, 0.0)
    """

    def __init__(self, *args):
        super(Direction, self).__init__(*args)
        Geometry.__init__(self)

    def __str__(self):
        return 'Direction({0}, {1}, {2})'.format(*self.xyz)

    def __repr__(self):
        return 'Direction({0}, {1}, {2})'.format(*self.xyz)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.ijk, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.ijk:
            yield elm

    def __len__(self):
        return 3

    def __getitem__(self, item):
        return self.ijk[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def i(self):
        """
        The direction i-component.

        :getter: Returns the i-component.
        :setter: Sets the i-component.
        :type: float
        """
        return self.X()

    @i.setter
    def i(self, i):
        self.SetX(i)

    @property
    def j(self):
        """
        The direction j-component.

        :getter: Returns the j-component.
        :setter: Sets the j-component.
        :type: float
        """
        return self.Y()

    @j.setter
    def j(self, j):
        self.SetY(j)

    @property
    def k(self):
        """
        The direction k-component.

        :getter: Returns the k-component.
        :setter: Sets the k-component.
        :type: float
        """
        return self.Z()

    @k.setter
    def k(self, k):
        self.SetZ(k)

    @property
    def ijk(self):
        """
        :return: The direction ijk-components.
        :rtype: numpy.ndarray
        """
        return array([self.i, self.j, self.k], dtype=float64)

    @property
    def xyz(self):
        """
        :return: The direction ijk-components (Same as ijk property,
            for compatibility only).
        :rtype: numpy.ndarray
        """
        return self.ijk

    @property
    def mag(self):
        """
        :return: Direction magnitude is always 1.
        :rtype: float
        """
        return 1.


class Vector(gp_Vec, Geometry):
    """
    Vector in 3-D space. Supports NumPy array methods.

    For more information see gp_Vec_.

    .. _gp_Vec: https://www.opencascade.com/doc/occt-7.1.0/refman/htmlclassgp___vec.html

    Usage:

    >>> from afem.geometry import Direction, Point, Vector
    >>> Vector(1., 2., 3.)
    Vector(1.0, 2.0, 3.0)
    >>> d = Direction(1., 0., 0.)
    >>> Vector(d)
    Vector(1.0, 0.0, 0.0)
    >>> p1 = Point()
    >>> p2 = Point(1., 2., 3.)
    >>> Vector(p1, p2)
    Vector(1.0, 2.0, 3.0)
    """

    def __init__(self, *args):
        super(Vector, self).__init__(*args)
        Geometry.__init__(self)

    def __str__(self):
        return 'Vector({0}, {1}, {2})'.format(*self.xyz)

    def __repr__(self):
        return 'Vector({0}, {1}, {2})'.format(*self.xyz)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self.xyz, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self.xyz:
            yield elm

    def __len__(self):
        return 3

    def __getitem__(self, item):
        return self.xyz[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def x(self):
        """
        The vector x-component.

        :getter: Returns the x-component.
        :setter: Sets the x-component.
        :type: float
        """
        return self.X()

    @x.setter
    def x(self, x):
        self.SetX(x)

    @property
    def y(self):
        """
        The vector y-component.

        :getter: Returns the y-component.
        :setter: Sets the y-component.
        :type: float
        """
        return self.Y()

    @y.setter
    def y(self, y):
        self.SetY(y)

    @property
    def z(self):
        """
        The vector z-component.

        :getter: Returns the z-component.
        :setter: Sets the z-component.
        :type: float
        """
        return self.Z()

    @z.setter
    def z(self, z):
        self.SetZ(z)

    @property
    def xyz(self):
        """
        :return: The vector xyz-components.
        :rtype: numpy.ndarray
        """
        return array([self.x, self.y, self.z], dtype=float64)

    @property
    def mag(self):
        """
        :return: Vector magnitude.
        :rtype: float
        """
        return self.Magnitude()

    @property
    def ijk(self):
        """
        :return: Normalized vector xyz-components.
        """
        return self.xyz / self.mag

    def reverse(self):
        """
        Reverse the direction of the vector.

        :return: None.
        """
        self.Reverse()

    def normalize(self):
        """
        Normalize the vector.

        :return: None.
        """
        self.Normalize()

    def scale(self, scale):
        """
        Scale the vector.

        :param float scale: Scaling value.

        :return: None.
        """
        self.Scale(scale)


class Axis1(gp_Ax1):
    """
    Axis in 3-D space. Definition incomplete.

    For more information see gp_Ax1_.

    .. _gp_Ax1: https://www.opencascade.com/doc/occt-7.1.0/refman/html/classgp___ax1.html

    Usage:

    >>> from afem.geometry import Axis1, Direction, Point
    >>> ax1 = Axis1()
    >>> p = Point()
    >>> d = Direction(1., 0., 0.)
    >>> ax1_ = Axis1(p, d)
    """

    def __init__(self, *args):
        super(Axis1, self).__init__(*args)


class Axis3(gp_Ax3):
    """
    Coordinate system in 3-D space. Definition incomplete.

    For more information see gp_Ax3_.

    .. _gp_Ax3: https://www.opencascade.com/doc/occt-7.1.0/refman/html/classgp___ax3.html

    Usage:

    >>> from afem.geometry import Axis3, Direction, Point
    >>> p = Point()
    >>> n = Direction(0., 0., 1.)
    >>> vx = Direction(1., 0., 0.)
    >>> ax3 = Axis3(p, n, vx)
    >>> ax3_ = Axis3(p, n)
    """

    def __init__(self, *args):
        super(Axis3, self).__init__(*args)


class Curve(Geometry):
    """
    Base class for curves.

    :param OCC.Geom.Handle_Geom_Curve the_handle: The curve handle.

    For more information see Geom_Curve_.

    .. _Geom_Curve: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___curve.html
    """

    def __init__(self, the_handle):
        super(Curve, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Curve
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_Curve
        """
        return self._object

    @property
    def adaptor(self):
        """
        :return: A curve adaptor.
        :rtype: OCC.GeomAdaptor.GeomAdaptor_Curve
        """
        return GeomAdaptor_Curve(self.handle)

    @property
    def u1(self):
        """
        :return: The first parameter.
        :rtype: float
        """
        return self.object.FirstParameter()

    @property
    def u2(self):
        """
        :return: The last parameter.
        :rtype: float
        """
        return self.object.LastParameter()

    @property
    def is_closed(self):
        """
        :return: *True* if curve is closed, *False* if not.
        :rtype: bool
        """
        return self.object.IsClosed()

    @property
    def is_periodic(self):
        """
        :return: *True* if curve is periodic, *False* if not.
        :rtype: bool
        """
        return self.object.IsPeriodic()

    @property
    def p1(self):
        """
        :return: The first point.
        :rtype: afem.geometry.entities.Point
        """
        return self.eval(self.u1)

    @property
    def p2(self):
        """
        :return: The last point.
        :rtype: afem.geometry.entities.Point
        """
        return self.eval(self.u2)

    @property
    def length(self):
        """
        :return: Curve length.
        :rtype: float
        """
        return self.arc_length(self.u1, self.u2)

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.Curve
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return Curve(h_crv)

    def local_to_global_param(self, *args):
        """
        Convert parameter(s) from local domain 0. <= u <= 1. to global domain
        a <= u <= b.

        :param float args: Local parameter(s).

        :return: Global parameter(s).
        :rtype: float or list[float]
        """
        return local_to_global_param(self.u1, self.u2, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0. <= u <= 1.

        :param float args: Global parameter(s).

        :return: Local parameter(s).
        :rtype: float or list[float]
        """
        return global_to_local_param(self.u1, self.u2, *args)

    def eval(self, u):
        """
        Evaluate a point on the curve.

        :param float u: Curve parameter.

        :return: Curve point.
        :rtype: afem.geometry.entities.Point
        """
        p = Point()
        self.object.D0(u, p)
        return p

    def deriv(self, u, d=1):
        """
        Evaluate a derivative on the curve.

        :param float u: Curve parameter.
        :param int d: Derivative to evaluate.

        :return: Curve derivative.
        :rtype: afem.geometry.entities.Vector
        """
        return Vector(self.object.DN(u, d).XYZ())

    def reverse(self):
        """
        Reverse curve direction.

        :return: None.
        """
        self.object.Reverse()

    def reversed_u(self, u):
        """
        Calculate the parameter on the reversed curve.

        :param float u: Curve parameter.

        :return: Reversed parameter.
        :rtype: float
        """
        return self.object.ReversedParameter(u)

    def arc_length(self, u1, u2, tol=1.0e-7):
        """
        Calculate the curve length between the parameters.

        :param float u1: First parameter.
        :param float u2: Last parameter.
        :param float tol: The tolerance.

        :return: Curve length.
        :rtype: float
        """
        if u1 > u2:
            u1, u2 = u2, u1
        return GCPnts_AbscissaPoint.Length(self.adaptor, u1, u2, tol)


class Line(Curve):
    """
    Infinite line.

    :param OCC.Geom.Handle_Geom_Line the_handle: The line handle.

    For more information see Geom_Line_.

    .. _Geom_Line: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___line.html
    """

    def __init__(self, the_handle):
        super(Line, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Line
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_Line
        """
        return self._object

    @staticmethod
    def downcast(crv):
        """
        Downcast the curve to this type.

        :param afem.geometry.entities.Curve crv: The curve.

        :return: A line.
        :rtype: afem.geometry.entities.Line

        :raise ValueError: If the curve cannot be downcast to this type.
        """
        h_crv = Handle_Geom_Line.DownCast(crv.handle)
        if h_crv.IsNull():
            raise ValueError('Could not downcast curve.')
        return Line(h_crv)

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.Line
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return Line(h_crv)


class Circle(Curve):
    """
    Circular curve.

    :param OCC.Geom.Handle_Geom_Circle the_handle: The circle handle.
    """

    def __init__(self, the_handle):
        super(Circle, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Circle
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_Circle
        """
        return self._object

    @property
    def radius(self):
        """
        :return: The radius.
        :rtype: float
        """
        return self.object.Radius()

    @staticmethod
    def downcast(crv):
        """
        Downcast the curve to this type.

        :param afem.geometry.entities.Curve crv: The curve.

        :return: A circle.
        :rtype: afem.geometry.entities.Circle

        :raise ValueError: If the curve cannot be downcast to this type.
        """
        h_crv = Handle_Geom_Circle.DownCast(crv.handle)
        if h_crv.IsNull():
            raise ValueError('Could not downcast curve.')
        return Circle(h_crv)

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.Circle
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return Circle(h_crv)

    def set_radius(self, r):
        """
        Set the radius.

        :param float r: The radius.

        :return: None.
        """
        self.object.SetRadius(r)


class Ellipse(Curve):
    """
    Elliptical curve.

    :param OCC.Geom.Handle_Geom_Ellipse the_handle: The ellipse handle.
    """

    def __init__(self, the_handle):
        super(Ellipse, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Ellipse
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_Ellipse
        """
        return self._object

    @property
    def major_radius(self):
        """
        :return: The major radius.
        :rtype: float
        """
        return self.object.MajorRadius()

    @property
    def minor_radius(self):
        """
        :return: The minor radius.
        :rtype: float
        """
        return self.object.MinorRadius()

    @staticmethod
    def downcast(crv):
        """
        Downcast the curve to this type.

        :param afem.geometry.entities.Curve crv: The curve.

        :return: An ellipse.
        :rtype: afem.geometry.entities.Ellipse

        :raise ValueError: If the curve cannot be downcast to this type.
        """
        h_crv = Handle_Geom_Ellipse.DownCast(crv.handle)
        if h_crv.IsNull():
            raise ValueError('Could not downcast curve.')
        return Ellipse(h_crv)

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.Ellipse
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return Ellipse(h_crv)

    def set_major_radius(self, r):
        """
        Set the major radius.

        :param float r: The major radius.

        :return: None.
        """
        self.object.SetMajorRadius(r)

    def set_minor_radius(self, r):
        """
        Set the minor radius.

        :param float r: The minor radius.

        :return: None.
        """
        self.object.SetMinorRadius(r)


class NurbsCurve(Curve):
    """
    NURBS curve in 3-D space.

    :param OCC.Geom.Handle_Geom_BSplineCurve the_handle: The curve handle.

    For more information see Geom_BSplineCurve_.

    .. _Geom_BSplineCurve: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___b_spline_curve.html
    """

    def __init__(self, the_handle):
        super(NurbsCurve, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_BSplineCurve
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_BSplineCurve
        """
        return self._object

    @property
    def p(self):
        """
        :return: Degree of curve.
        :rtype: int
        """
        return self.object.Degree()

    @property
    def n(self):
        """
        :return: Number of control points.
        :rtype: int
        """
        return self.object.NbPoles()

    @property
    def knots(self):
        """
        :return: Knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.object.NbKnots())
        self.object.Knots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def mult(self):
        """
        :return: Multiplicity of knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.object.NbKnots())
        self.object.Multiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        """
        :return: Knot sequence.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.object.NbPoles() +
                                             self.object.Degree() + 1)
        self.object.KnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def cp(self):
        """
        :return: Control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColgp_Array1OfPnt(1, self.object.NbPoles())
        self.object.Poles(tcol_array)
        return to_np_from_tcolgp_array1_pnt(tcol_array)

    @property
    def w(self):
        """
        :return: Weights of control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.object.NbPoles())
        self.object.Weights(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def cpw(self):
        """
        :return: Homogeneous control points.
        :rtype: numpy.ndarray
        """
        return homogenize_array1d(self.cp, self.w)

    @staticmethod
    def downcast(crv):
        """
        Downcast the curve to this type.

        :param afem.geometry.entities.Curve crv: The curve.

        :return: A NURBS curve.
        :rtype: afem.geometry.entities.NurbsCurve

        :raise ValueError: If the curve cannot be downcast to this type.
        """
        h_crv = Handle_Geom_BSplineCurve.DownCast(crv.handle)
        if h_crv.IsNull():
            raise ValueError('Could not downcast curve.')
        return NurbsCurve(h_crv)

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.NurbsCurve
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return NurbsCurve(h_crv)

    def set_domain(self, u1=0., u2=1.):
        """
        Reparameterize the knot vector between *u1* and *u2*.

        :param float u1: First parameter.
        :param float u2: Last parameter.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if u1 > u2:
            return False
        tcol_knots = TColStd_Array1OfReal(1, self.object.NbKnots())
        self.object.Knots(tcol_knots)
        reparameterize_knots(u1, u2, tcol_knots)
        self.object.SetKnots(tcol_knots)
        return True

    def segment(self, u1, u2):
        """
        Segment the curve between parameters.

        :param float u1: First parameter.
        :param float u2: Last parameter.

        :return: *True* if segmented, *False* if not.
        :rtype: bool
        """
        if u1 > u2:
            return False
        self.object.Segment(u1, u2)
        return True

    def set_cp(self, i, cp, weight=None):
        """
        Modify the curve by setting the specified control point and weight.

        :param int i: The point index (1 <= *i* <= *n*).
        :param afem.geometry.entities.Point cp: The point.
        :param weight: The weight.
        :type weight: float or None

        :return: None.
        """
        if weight is None:
            self.object.SetPole(i, cp)
        else:
            self.object.SetPole(i, cp, weight)


class TrimmedCurve(Curve):
    """
    Trimmed curve. This defines a basis curve limited by two parameter values.

    :param OCC.Geom.Handle_Geom_TrimmedCurve the_handle: The handle.
    """

    def __init__(self, the_handle):
        super(TrimmedCurve, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_TrimmedCurve
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_TrimmedCurve
        """
        return self._object

    @property
    def basis_curve(self):
        """
        :return: The basis curve.
        :rtype: afem.geometry.entities.Curve
        """
        return Curve(self.object.BasisCurve())

    @staticmethod
    def downcast(crv):
        """
        Downcast the curve to this type.

        :param afem.geometry.entities.Curve crv: The curve.

        :return: A trimmed curve.
        :rtype: afem.geometry.entities.TrimmedCurve

        :raise ValueError: If the curve cannot be downcast to this type.
        """
        h_crv = Handle_Geom_TrimmedCurve.DownCast(crv.handle)
        if h_crv.IsNull():
            raise ValueError('Could not downcast curve.')
        return TrimmedCurve(h_crv)

    def copy(self):
        """
        Return a new copy of the curve.

        :return: New curve.
        :rtype: afem.geometry.entities.Circle
        """
        h_geom = self.object.Copy()
        h_crv = self.handle.DownCast(h_geom)
        return TrimmedCurve(h_crv)

    def set_trim(self, u1, u2, sense=True, adjust_periodic=True):
        """
        Set the trimming parameters on the basis curve.

        :param float u1: The first parameter.
        :param float u2: The last parameter.
        :param bool sense: If the basis curve is periodic, the trimmed curve
            will have the same orientation as the basis curve if ``True`` or
            opposite if ``False``.
        :param bool adjust_periodic: If the basis curve is periodic, the bounds
            of the trimmed curve may be different from *u1* and *u2* if
            ``True``.

        :return: None.

        :raise RuntimeError: If *u1* == *u2*.
        :raise RuntimeError: If *u1* or *u2* is outside the bounds of the basis
            curve.
        """
        self.object.SetTrim(u1, u2, sense, adjust_periodic)


class Surface(Geometry):
    """
    Base class for surfaces.

    :param OCC.Geom.Handle_Geom_Surface the_handle: The surface handle.

    For more information see Geom_Surface_.

    .. _Geom_Surface: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___surface.html
    """

    def __init__(self, the_handle):
        super(Surface, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Surface
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_Surface
        """
        return self._object

    @property
    def adaptor(self):
        """
        :return: A surface adaptor.
        :rtype: OCC.GeomAdaptor.GeomAdaptor_Surface
        """
        return GeomAdaptor_Surface(self.handle)

    @property
    def u1(self):
        """
        :return: The first parameter in u-direction.
        :rtype: float
        """
        return self.object.Bounds()[0]

    @property
    def u2(self):
        """
        :return: The last parameter in u-direction.
        :rtype: float
        """
        return self.object.Bounds()[1]

    @property
    def v1(self):
        """
        :return: The first parameter in v-direction.
        :rtype: float
        """
        return self.object.Bounds()[2]

    @property
    def v2(self):
        """
        :return: The last parameter in v-direction.
        :rtype: float
        """
        return self.object.Bounds()[3]

    @property
    def area(self):
        """
        :return: The surface area.
        :rtype: float
        """
        return self.surface_area(self.u1, self.v1, self.u2, self.v1)

    def copy(self):
        """
        Return a new copy of the surface.

        :return: New surface.
        :rtype: afem.geometry.entities.Surface
        """
        h_geom = self.object.Copy()
        h_srf = self.handle.DownCast(h_geom)
        return Surface(h_srf)

    def eval(self, u=0., v=0.):
        """
        Evaluate a point on the surface.

        :param float u: Surface u-parameter.
        :param float v: Surface v-parameter.

        :return: Surface point.
        :rtype: afem.geometry.entities.Point
        """
        p = Point()
        self.object.D0(u, v, p)
        return p

    def deriv(self, u, v, nu, nv):
        """
        Evaluate a derivative on the surface.

        :param float u: Surface u-parameter.
        :param float v: Surface v-parameter.
        :param int nu: Derivative in u-direction.
        :param int nv: Derivative in v-direction.

        :return: Surface derivative.
        :rtype: afem.geometry.entities.Vector
        """
        return Vector(self.object.DN(u, v, nu, nv).XYZ())

    def norm(self, u, v):
        """
        Evaluate a normal on the surface.

        :param float u: Surface u-parameter.
        :param float v: Surface v-parameter.

        :return: Surface normal.
        :rtype: afem.geometry.entities.Vector
        """
        du = self.deriv(u, v, 1, 0)
        dv = self.deriv(u, v, 0, 1)
        return Vector(du.Crossed(dv).XYZ())

    def surface_area(self, u1, v1, u2, v2, tol=1.0e-7):
        """
        Calculate the surface area between the parameter.

        :param float u1:
        :param float v1:
        :param float u2:
        :param float v2:
        :param float tol: The tolerance.

        :return: The area.
        :rtype: float
        """
        f = BRepBuilderAPI_MakeFace(self.handle, u1, u2, v1, v2, tol).Face()
        sprops = GProp_GProps()
        brepgprop_SurfaceProperties(f, sprops, tol)
        return sprops.Mass()

    def u_iso(self, u):
        """
        Get a iso-parametric curve at a constant u-parameter.

        :param float u: The u-parameter.

        :return: The curve.
        :rtype: afem.geometry.entities.Curve
        """
        h_crv = self.object.UIso(u)
        return Curve(h_crv)

    def v_iso(self, v):
        """
        Get a iso-parametric curve at a constant v-parameter.

        :param float v: The v-parameter.

        :return: The curve.
        :rtype: afem.geometry.entities.Curve
        """
        h_crv = self.object.VIso(v)
        return Curve(h_crv)


class Plane(Surface):
    """
    Infinite plane.

    :param OCC.Geom.Handle_Geom_Plane the_handle: The plane handle.

    For more information see Geom_Plane_.

    .. _Geom_Plane: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___plane.html

    """

    def __init__(self, the_handle):
        super(Plane, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Plane
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_Plane
        """
        return self._object

    @staticmethod
    def downcast(srf):
        """
        Downcast the surface to this type.

        :param afem.geometry.entities.Surface srf: The surface.

        :return: A plane.
        :rtype: afem.geometry.entities.Plane

        :raise ValueError: If the surface cannot be downcast to this type.
        """
        h_srf = Handle_Geom_Plane.DownCast(srf.handle)
        if h_srf.IsNull():
            raise ValueError('Could not downcast surface.')
        return Plane(h_srf)

    def copy(self):
        """
        Return a new copy of the plane.

        :return: New plane.
        :rtype: afem.geometry.entities.Plane
        """
        h_geom = self.object.Copy()
        h_srf = self.handle.DownCast(h_geom)
        return Plane(h_srf)

    def distance(self, pnt):
        """
        Compute the distance between a point and this plane.

        :param point_like pnt: The point.

        :return: The distance.
        :rtype: float

        :raises TypeError: If *pnt* cannot be converted to a point.
        """
        if isinstance(pnt, gp_Pnt):
            pnt = Point(pnt.XYZ())
        elif isinstance(pnt, (tuple, list, ndarray)):
            pnt = Point(pnt[0], pnt[1], pnt[2])
        else:
            msg = 'Invalid point type.'
            raise TypeError(msg)

        return self.object.Pln().Distance(pnt)


class NurbsSurface(Surface):
    """
    NURBS surface in 3-D space.

    :param OCC.Geom.Handle_Geom_BSplineSurface the_handle: The surface handle.

    For more information see Geom_BSplineSurface_.

    .. _Geom_BSplineSurface: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___b_spline_surface.html

    """

    def __init__(self, the_handle):
        super(NurbsSurface, self).__init__(the_handle)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_BSplineSurface
        """
        return self._handle

    @property
    def object(self):
        """
        :return: The underlying object.
        :rtype: OCC.Geom.Geom_BSplineSurface
        """
        return self._object

    @property
    def p(self):
        """
        :return: Degree in u-direction.
        :rtype: int
        """
        return self.object.UDegree()

    @property
    def q(self):
        """
        :return: Degree in v-direction.
        :rtype: int
        """
        return self.object.VDegree()

    @property
    def n(self):
        """
        :return: Number of control points in u-direction.
        :rtype: int
        """
        return self.object.NbUPoles()

    @property
    def m(self):
        """
        :return: Number of control points in v-direction.
        :rtype: int
        """
        return self.object.NbVPoles()

    @property
    def uknots(self):
        """
        :return: Knot vector in u-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.object.NbUKnots())
        self.object.UKnots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def umult(self):
        """
        :return: Multiplicity of knot vector in u-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.object.NbUKnots())
        self.object.UMultiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        """
        :return: Knot sequence in u-direction.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.object.NbUPoles() +
                                             self.object.UDegree() + 1)
        self.object.UKnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def vknots(self):
        """
        :return: Knot vector in v-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.object.NbVKnots())
        self.object.VKnots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def vmult(self):
        """
        :return: Multiplicity of knot vector in v-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.object.NbVKnots())
        self.object.VMultiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def vk(self):
        """
        :return: Knot sequence in v-direction.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.object.NbVPoles() +
                                             self.object.VDegree() + 1)
        self.object.VKnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def cp(self):
        """
        :return: Control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColgp_Array2OfPnt(1, self.object.NbUPoles(),
                                        1, self.object.NbVPoles())
        self.object.Poles(tcol_array)
        return to_np_from_tcolgp_array2_pnt(tcol_array)

    @property
    def w(self):
        """
        :return: Weights of control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array2OfReal(1, self.object.NbUPoles(),
                                          1, self.object.NbVPoles())
        self.object.Weights(tcol_array)
        return to_np_from_tcolstd_array2_real(tcol_array)

    @property
    def cpw(self):
        """
        :return: Homogeneous control points.
        :rtype: numpy.ndarray
        """
        return homogenize_array2d(self.cp, self.w)

    @staticmethod
    def downcast(srf):
        """
        Downcast the surface to this type.

        :param afem.geometry.entities.Surface srf: The surface.

        :return: A surface.
        :rtype: afem.geometry.entities.NurbsSurface

        :raise ValueError: If the surface cannot be downcast to this type.
        """
        h_srf = Handle_Geom_BSplineSurface.DownCast(srf.handle)
        if h_srf.IsNull():
            raise ValueError('Could not downcast surface.')
        return NurbsSurface(h_srf)

    def copy(self):
        """
        Return a new copy of the surface.

        :return: New surface.
        :rtype: afem.geometry.entities.NurbsSurface
        """
        h_geom = self.object.Copy()
        h_srf = self.handle.DownCast(h_geom)
        return NurbsSurface(h_srf)

    def set_udomain(self, u1=0., u2=1.):
        """
        Reparameterize the knot vector between *u1* and *u2*.

        :param float u1: First parameter.
        :param float u2: Last parameter.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if u1 > u2:
            return False
        tcol_knots = TColStd_Array1OfReal(1, self.object.NbUKnots())
        self.object.UKnots(tcol_knots)
        reparameterize_knots(u1, u2, tcol_knots)
        self.object.SetUKnots(tcol_knots)
        return True

    def set_vdomain(self, v1=0., v2=1.):
        """
        Reparameterize the knot vector between *v1* and *v2*.

        :param float v1: First parameter.
        :param float v2: Last parameter.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if v1 > v2:
            return False
        tcol_knots = TColStd_Array1OfReal(1, self.object.NbVKnots())
        self.object.VKnots(tcol_knots)
        reparameterize_knots(v1, v2, tcol_knots)
        self.object.SetVKnots(tcol_knots)
        return True

    def local_to_global_param(self, d, *args):
        """
        Convert parameter(s) from local domain 0. <= u,v <= 1. to global domain
        a <= u,v <= b.

        :param str d: Parameter direction ('u' or 'v').
        :param float args: Local parameter(s).

        :return: Global parameter(s).
        :rtype: float or list[float]
        """
        if d.lower() in ['u']:
            return local_to_global_param(self.u1, self.u2, *args)
        else:
            return local_to_global_param(self.v1, self.v2, *args)

    def global_to_local_param(self, d, *args):
        """
        Convert surface parameter(s) from global domain a <= u,v <= b to local
        domain 0. <= u,v <= 1.

        :param str d: Parameter direction ('u' or 'v').
        :param float args: Global parameter(s).

        :return: Local parameter(s).
        :rtype: float or list[float]
        """
        if d.lower() in ['u']:
            return global_to_local_param(self.u1, self.u2, *args)
        else:
            return global_to_local_param(self.v1, self.v2, *args)

    def segment(self, u1, u2, v1, v2):
        """
        Segment the surface between parameters.

        :param float u1: First parameter in u-direction.
        :param float u2: Last parameter in u-direction.
        :param float v1: First parameter in v-direction.
        :param float v2: Last parameter in v-direction.

        :return: *True* if segmented, *False* if not.
        :rtype: bool
        """
        if u1 > u2 or v1 > v2:
            return False
        self.object.CheckAndSegment(u1, u2, v1, v2)
        return True

    def locate_u(self, u, tol2d=1.0e-9, with_knot_repetition=False):
        """
        Locate the u-value in the knot sequence.

        :param float u: The parameter.
        :param int tol2d: The parametric tolerance. Used to determine if *u* is
            at an existing knot.
        :param bool with_knot_repetition: Considers location of knot value
            with repetition of multiple knot value if ``True``.

        :return: Bounding knot locations (i1, i2).
        :rtype: tuple(float)
        """
        return self.object.LocateU(u, tol2d, with_knot_repetition)

    def locate_v(self, v, tol2d=1.0e-9, with_knot_repetition=False):
        """
        Locate the v-value in the knot sequence.

        :param float v: The parameter.
        :param int tol2d: The parametric tolerance. Used to determine if *v* is
            at an existing knot.
        :param bool with_knot_repetition: Considers location of knot value
            with repetition of multiple knot value if ``True``.

        :return: Bounding knot locations (i1, i2).
        :rtype: tuple(float)
        """
        return self.object.LocateV(v, tol2d, with_knot_repetition)

    def insert_uknot(self, u, m=1, tol2d=1.0e-9):
        """
        Insert a u-value in the knot sequence.

        :param float u: The knot value.
        :param int m: If the knot already exists, the multiplicity of the
            knot is increased if the previous multiplicity is lower than
            *m*. Otherwise it does nothing.
        :param float tol2d: Parametric tolerance for comparing knot values.

        :return: None.
        """
        self.object.InsertUKnot(u, m, tol2d)

    def insert_vknot(self, v, m=1, tol2d=1.0e-9):
        """
        Insert a v-value in the knot sequence.

        :param float v: The knot value.
        :param int m: If the knot already exists, the multiplicity of the
            knot is increased if the previous multiplicity is lower than
            *m*. Otherwise it does nothing.
        :param float tol2d: Parametric tolerance for comparing knot values.

        :return: None.
        """
        self.object.InsertVKnot(v, m, tol2d)

    def set_cp(self, i, j, cp, weight=None):
        """
        Modify the surface by setting the specified control point and weight.

        :param int i: The point index in u-direction (1 <= *i* <= *n*).
        :param int j: The point index in v-direction (1 <= *j* <= *m*).
        :param afem.geometry.entities.Point cp: The point.
        :param weight: The weight.
        :type weight: float or None

        :return: None.
        """
        if weight is None:
            self.object.SetPole(i, j, cp)
        else:
            self.object.SetPole(i, j, cp, weight)

    def set_cp_row(self, u_index, cp, weights=None):
        """
        Modify the surface by setting the specified row of control points
        and weights.

        :param int u_index: The row index (1 <= *u_index* <= *n*).
        :param list[afem.geometry.entities.Point] cp: The points.
        :param weights: The weights.
        :type weights: list[float] or None

        :return: None.
        """
        tcol_gp = to_tcolgp_array1_pnt(cp)
        if weights is None:
            self.object.SetPoleRow(u_index, tcol_gp)
        else:
            tcol_w = to_tcolstd_array1_real(weights)
            self.object.SetPoleRow(u_index, tcol_gp, tcol_w)

    def set_cp_col(self, v_index, cp, weights=None):
        """
        Modify the surface by setting the specified column of control points
        and weights.

        :param int v_index: The column index (1 <= *v_index* <= *m*).
        :param list[afem.geometry.entities.Point] cp: The points.
        :param weights: The weights.
        :type weights: list[float] or None

        :return: None.
        """
        tcol_gp = to_tcolgp_array1_pnt(cp)
        if weights is None:
            self.object.SetPoleCol(v_index, tcol_gp)
        else:
            tcol_w = to_tcolstd_array1_real(weights)
            self.object.SetPoleCol(v_index, tcol_gp, tcol_w)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
