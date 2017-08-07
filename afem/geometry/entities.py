from OCC.Geom import (Geom_BSplineCurve, Geom_BSplineSurface, Geom_Curve,
                      Geom_Line, Geom_Plane, Geom_Surface)
from OCC.Geom2d import Geom2d_BSplineCurve
from OCC.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface
from OCC.TColStd import (TColStd_Array1OfInteger, TColStd_Array1OfReal,
                         TColStd_Array2OfReal)
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.gp import gp_Ax1, gp_Ax3, gp_Dir, gp_Pnt, gp_Pnt2d, gp_Vec, gp_XYZ
from numpy import add, array, float64, ndarray, subtract

from afem.config import Settings
from afem.geometry.methods.calculate import curve_length
from afem.geometry.methods.parameterize import reparameterize_knots
from afem.geometry.utils import (global_to_local_param, homogenize_array1d,
                                 homogenize_array2d, local_to_global_param)
from afem.graphics.viewer import ViewableItem
from afem.occ.utils import (to_np_from_tcolgp_array1_pnt,
                            to_np_from_tcolgp_array2_pnt,
                            to_np_from_tcolstd_array1_integer,
                            to_np_from_tcolstd_array1_real,
                            to_np_from_tcolstd_array2_real)
from afem.utils.check import is_array_like

__all__ = ["Geometry", "Point", "Point2D", "Direction", "Vector", "Axis1",
           "Axis3", "Curve", "Line", "NurbsCurve", "NurbsCurve2D", "Surface",
           "Plane", "NurbsSurface"]


class Geometry(ViewableItem):
    """
    Base class for geometry.
    """

    def __init__(self):
        super(Geometry, self).__init__()


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

    def copy(self):
        """
        Return a new copy of the point.

        :return: New point.
        :rtype: afem.geometry.entities.Point
        """
        return Point(*self.xyz)

    def set_xyz(self, xyz):
        """
        Set point coordinates.

        :param array_like xyz: Point coordinates.

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

    def is_equal(self, other, tol=None):
        """
        Check for coincident points.

        :param point_like other: The other point.
        :param float tol: Tolerance for coincidence.

        :return: *True* if coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.gtol
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
            return Vector(v[0], v[1], v[2])
        elif isinstance(v, Direction):
            v = Vector(v)
        elif isinstance(v, gp_Vec):
            v = Vector(v.XYZ())
        else:
            msg = 'Invalid vector type.'
            raise TypeError(msg)

        self.Translate(v)
        return True


class Point2D(gp_Pnt2d, Geometry):
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
        Geometry.__init__(self)

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

    def copy(self):
        """
        Return a new copy of the point.

        :return: New point.
        :rtype: afem.geometry.entities.Point2D
        """
        return Point2D(*self.xy)

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

    def is_equal(self, other, tol=None):
        """
        Check for coincident points.

        :param point_like other: The other point.
        :param float tol: Tolerance for coincidence.

        :return: *True* if coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.gtol
        if isinstance(other, gp_Pnt2d):
            return self.IsEqual(other, tol)
        if is_array_like(other) and len(other) == 2:
            other = Point2D(*other)
            return self.IsEqual(other, tol)
        return False


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


class Curve(Geom_Curve, Geometry):
    """
    Base class for curves.

    For more information see Geom_Curve_.

    .. _Geom_Curve: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___curve.html
    """

    def __init__(self):
        Geometry.__init__(self)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Curve
        """
        return self.GetHandle()

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
        return self.FirstParameter()

    @property
    def u2(self):
        """
        :return: The last parameter.
        :rtype: float
        """
        return self.LastParameter()

    @property
    def is_closed(self):
        """
        :return: *True* if curve is closed, *False* if not.
        :rtype: bool
        """
        return self.IsClosed()

    @property
    def is_periodic(self):
        """
        :return: *True* if curve is periodic, *False* if not.
        :rtype: bool
        """
        return self.IsPeriodic()

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

    def eval(self, u):
        """
        Evaluate a point on the curve.

        :param float u: Curve parameter.

        :return: Curve point.
        :rtype: afem.geometry.entities.Point
        """
        p = Point()
        self.D0(u, p)
        return p

    def deriv(self, u, d=1):
        """
        Evaluate a derivative on the curve.

        :param float u: Curve parameter.
        :param int d: Derivative to evaluate.

        :return: Curve derivative.
        :rtype: afem.geometry.entities.Vector
        """
        return Vector(self.DN(u, d).XYZ())

    def reverse(self):
        """
        Reverse curve direction.

        :return: None.
        """
        self.Reverse()

    def reversed_u(self, u):
        """
        Calculate the parameter on the reversed curve.

        :param float u: Curve parameter.

        :return: Reversed parameter.
        :rtype: float
        """
        return self.ReversedParameter(u)

    def arc_length(self, u1, u2):
        """
        Calculate the curve length between the parameters.

        :param float u1: First parameter.
        :param float u2: Last parameter.

        :return: Curve length.
        :rtype: float
        """
        return curve_length(self, u1, u2)


class Line(Geom_Line, Curve):
    """
    Infinite line.

    For more information see Geom_Line_.

    .. _Geom_Line: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___line.html

    Usage:

    >>> from afem.geometry import Direction, Line, Point
    >>> p = Point()
    >>> d = Direction(1., 0., 0.)
    >>> line = Line(p, d)
    """

    def __init__(self, *args):
        super(Line, self).__init__(*args)
        Curve.__init__(self)


class NurbsCurve(Geom_BSplineCurve, Curve):
    """
    NURBS curve in 3-D space.

    For more information see Geom_BSplineCurve_.

    .. _Geom_BSplineCurve: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___b_spline_curve.html
    """

    def __init__(self, *args):
        super(NurbsCurve, self).__init__(*args)
        Curve.__init__(self)

    @property
    def p(self):
        """
        :return: Degree of curve.
        :rtype: int
        """
        return self.Degree()

    @property
    def n(self):
        """
        :return: Number of control points - 1.
        :rtype: int
        """
        return self.NbPoles() - 1

    @property
    def knots(self):
        """
        :return: Knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.NbKnots())
        self.Knots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def mult(self):
        """
        :return: Multiplicity of knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.NbKnots())
        self.Multiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        """
        :return: Knot sequence.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbPoles() +
                                             self.Degree() + 1)
        self.KnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def cp(self):
        """
        :return: Control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColgp_Array1OfPnt(1, self.NbPoles())
        self.Poles(tcol_array)
        return to_np_from_tcolgp_array1_pnt(tcol_array)

    @property
    def w(self):
        """
        :return: Weights of control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.NbPoles())
        self.Weights(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def cpw(self):
        """
        :return: Homogeneous control points.
        :rtype: numpy.ndarray
        """
        return homogenize_array1d(self.cp, self.w)

    @property
    def length(self):
        """
        :return: Curve length.
        :rtype: float
        """
        return self.arc_length(self.u1, self.u2)

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
        tcol_knots = TColStd_Array1OfReal(1, self.NbKnots())
        self.Knots(tcol_knots)
        reparameterize_knots(u1, u2, tcol_knots)
        self.SetKnots(tcol_knots)
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
        self.Segment(u1, u2)
        return True


class NurbsCurve2D(Geom2d_BSplineCurve, Geometry):
    """
    NURBS curve in 2-D space.

    For more information see Geom2d_BSplineCurve_.

    .. _Geom2d_BSplineCurve: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom2d___b_spline_curve.html
    """

    def __init__(self, *args):
        super(NurbsCurve2D, self).__init__(*args)
        Geometry.__init__(self)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom2d.Handle_Geom2d_BSplineCurve
        """
        return self.GetHandle()

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
        return self.FirstParameter()

    @property
    def u2(self):
        """
        :return: The last parameter.
        :rtype: float
        """
        return self.LastParameter()

    @property
    def p(self):
        """
        :return: Degree of curve.
        :rtype: int
        """
        return self.Degree()

    @property
    def n(self):
        """
        :return: Number of control points - 1.
        :rtype: int
        """
        return self.NbPoles() - 1

    @property
    def knots(self):
        """
        :return: Knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.NbKnots())
        self.Knots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def mult(self):
        """
        :return: Multiplicity of knot vector.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.NbKnots())
        self.Multiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        """
        :return: Knot sequence.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbPoles() +
                                             self.Degree() + 1)
        self.KnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def is_closed(self):
        """
        :return: *True* if curve is closed, *False* if not.
        :rtype: bool
        """
        return self.IsClosed()

    @property
    def is_periodic(self):
        """
        :return: *True* if curve is periodic, *False* if not.
        :rtype: bool
        """
        return self.IsPeriodic()

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
        tcol_knots = TColStd_Array1OfReal(1, self.NbKnots())
        self.Knots(tcol_knots)
        reparameterize_knots(u1, u2, tcol_knots)
        self.SetKnots(tcol_knots)
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
        self.D0(u, p)
        return p

    def reverse(self):
        """
        Reverse curve direction.

        :return: None.
        """
        self.Reverse()

    def reversed_u(self, u):
        """
        Calculate the parameter on the reversed curve.

        :param float u: Curve parameter.

        :return: Reversed parameter.
        :rtype: float
        """
        return self.ReversedParameter(u)

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
        self.Segment(u1, u2)
        return True


class Surface(Geom_Surface, Geometry):
    """
    Base class for surfaces.

    For more information see Geom_Surface_.

    .. _Geom_Surface: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___surface.html
    """

    def __init__(self):
        Geometry.__init__(self)

    @property
    def handle(self):
        """
        :return: The smart pointer.
        :rtype: OCC.Geom.Handle_Geom_Surface
        """
        return self.GetHandle()

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
        return self.Bounds()[0]

    @property
    def u2(self):
        """
        :return: The last parameter in u-direction.
        :rtype: float
        """
        return self.Bounds()[1]

    @property
    def v1(self):
        """
        :return: The first parameter in v-direction.
        :rtype: float
        """
        return self.Bounds()[2]

    @property
    def v2(self):
        """
        :return: The last parameter in v-direction.
        :rtype: float
        """
        return self.Bounds()[3]

    def eval(self, u=0., v=0.):
        """
        Evaluate a point on the surface.

        :param float u: Surface u-parameter.
        :param float v: Surface v-parameter.

        :return: Surface point.
        :rtype: afem.geometry.entities.Point
        """
        p = Point()
        self.D0(u, v, p)
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
        return Vector(self.DN(u, v, nu, nv).XYZ())

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


class Plane(Geom_Plane, Surface):
    """
    Infinite plane.

    For more information see Geom_Plane_.

    .. _Geom_Plane: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___plane.html

    """

    def __init__(self, *args):
        super(Plane, self).__init__(*args)
        Surface.__init__(self)


class NurbsSurface(Geom_BSplineSurface, Surface):
    """
    NURBS surface in 3-D space.

    For more information see Geom_BSplineSurface_.

    .. _Geom_BSplineSurface: https://www.opencascade.com/doc/occt-7.1.0/refman/html/class_geom___b_spline_surface.html

    """

    def __init__(self, *args):
        super(NurbsSurface, self).__init__(*args)
        Surface.__init__(self)

    @property
    def p(self):
        """
        :return: Degree in u-direction.
        :rtype: int
        """
        return self.UDegree()

    @property
    def q(self):
        """
        :return: Degree in v-direction.
        :rtype: int
        """
        return self.VDegree()

    @property
    def n(self):
        """
        :return: Number of control points - 1 in u-direction.
        :rtype: int
        """
        return self.NbUPoles() - 1

    @property
    def m(self):
        """
        :return: Number of control points - 1 in v-direction.
        :rtype: int
        """
        return self.NbVPoles() - 1

    @property
    def uknots(self):
        """
        :return: Knot vector in u-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.NbUKnots())
        self.UKnots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def umult(self):
        """
        :return: Multiplicity of knot vector in u-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.NbUKnots())
        self.UMultiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        """
        :return: Knot sequence in u-direction.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbUPoles() +
                                             self.UDegree() + 1)
        self.UKnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def vknots(self):
        """
        :return: Knot vector in v-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfReal(1, self.NbVKnots())
        self.VKnots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def vmult(self):
        """
        :return: Multiplicity of knot vector in v-direction.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array1OfInteger(1, self.NbVKnots())
        self.VMultiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def vk(self):
        """
        :return: Knot sequence in v-direction.
        :rtype: numpy.ndarray
        """
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbVPoles() +
                                             self.VDegree() + 1)
        self.VKnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def cp(self):
        """
        :return: Control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColgp_Array2OfPnt(1, self.NbUPoles(), 1, self.NbVPoles())
        self.Poles(tcol_array)
        return to_np_from_tcolgp_array2_pnt(tcol_array)

    @property
    def w(self):
        """
        :return: Weights of control points.
        :rtype: numpy.ndarray
        """
        tcol_array = TColStd_Array2OfReal(1, self.NbUPoles(),
                                          1, self.NbVPoles())
        self.Weights(tcol_array)
        return to_np_from_tcolstd_array2_real(tcol_array)

    @property
    def cpw(self):
        """
        :return: Homogeneous control points.
        :rtype: numpy.ndarray
        """
        return homogenize_array2d(self.cp, self.w)

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
        tcol_knots = TColStd_Array1OfReal(1, self.NbUKnots())
        self.UKnots(tcol_knots)
        reparameterize_knots(u1, u2, tcol_knots)
        self.SetUKnots(tcol_knots)
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
        tcol_knots = TColStd_Array1OfReal(1, self.NbVKnots())
        self.VKnots(tcol_knots)
        reparameterize_knots(v1, v2, tcol_knots)
        self.SetVKnots(tcol_knots)
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
        self.CheckAndSegment(u1, u2, v1, v2)
        return True


if __name__ == "__main__":
    import doctest

    doctest.testmod()
