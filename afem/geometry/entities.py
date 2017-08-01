from OCC.BRepBndLib import brepbndlib_Add
from OCC.Bnd import Bnd_Box
from OCC.Geom import (Geom_BSplineCurve, Geom_BSplineSurface, Geom_Curve,
                      Geom_Line, Geom_Plane, Geom_Surface)
from OCC.Geom2d import Geom2d_BSplineCurve
from OCC.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface
from OCC.TColStd import (TColStd_Array1OfInteger, TColStd_Array1OfReal,
                         TColStd_Array2OfReal)
from OCC.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt
from OCC.gp import gp_Ax1, gp_Ax3, gp_Dir, gp_Pnt, gp_Pnt2d, gp_Vec, gp_XYZ
from numpy import add, array, float64, subtract
from numpy.linalg import norm

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
           "Plane", "NurbsSurface", "BBox"]


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
    def xyz(self):
        """
        :return: The point location.
        :rtype: numpy.ndarray
        """
        return array([self.X(), self.Y(), self.Z()], dtype=float64)

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

        :param afem.geometry.entities.Vector v: The translation vector.

        :return: *True* if translated, *False* if not.
        :rtype: bool
        """
        self.Translate(v)
        return True


class Point2D(gp_Pnt2d, Geometry):
    """
    2-D point.
    """

    def __init__(self, *args):
        super(Point2D, self).__init__(*args)
        Geometry.__init__(self)

    def __str__(self):
        return 'Point 2-D = ({0}, {1})'.format(*self.xy)

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
        return array([self.X(), self.Y()], dtype=float64)

    @property
    def x(self):
        return self.X()

    @x.setter
    def x(self, x):
        self.SetX(x)

    @property
    def y(self):
        return self.Y()

    @y.setter
    def y(self, y):
        self.SetY(y)

    def copy(self):
        """
        Return a new copy of the point.

        :return:
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

        :param other: The other point.
        :type other: point_like

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

        :param other: The other point.
        :type other: point_like
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
    3-D unit vector.
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
        return self.X()

    @property
    def j(self):
        return self.Y()

    @property
    def k(self):
        return self.Z()

    @property
    def ijk(self):
        return array([self.i, self.j, self.k], dtype=float64)

    @property
    def xyz(self):
        return self.ijk

    @property
    def mag(self):
        return 1.


class Vector(gp_Vec, Geometry):
    """
    3-D vector.
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
        return self.X()

    @property
    def y(self):
        return self.Y()

    @property
    def z(self):
        return self.Z()

    @property
    def xyz(self):
        return array([self.x, self.y, self.z], dtype=float64)

    @property
    def mag(self):
        return norm(self.xyz)

    @property
    def ijk(self):
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

        :return:
        """
        self.Normalize()

    def scale(self, scale):
        """
        Scale the vector.

        :param scale:

        :return:
        """
        self.Scale(scale)


class Axis1(gp_Ax1):
    """
    Coordinate system.
    """

    def __init__(self, *args):
        super(Axis1, self).__init__(*args)


class Axis3(gp_Ax3):
    """
    Coordinate system.
    """

    def __init__(self, *args):
        super(Axis3, self).__init__(*args)


class Curve(Geom_Curve, Geometry):
    """
    Base class for curves.
    """

    def __init__(self):
        Geometry.__init__(self)

    @property
    def handle(self):
        return self.GetHandle()

    @property
    def adaptor(self):
        return GeomAdaptor_Curve(self.GetHandle())

    @property
    def u1(self):
        return self.FirstParameter()

    @property
    def u2(self):
        return self.LastParameter()

    @property
    def is_closed(self):
        return self.IsClosed()

    @property
    def is_periodic(self):
        return self.IsPeriodic()

    def eval(self, u):
        """
        Evaluate line at parameter.

        :param float u: Line parameter.

        :return: Point on line.
        :rtype: :class:`.Point`
        """
        p = Point()
        self.D0(u, p)
        return p

    def deriv(self, u, d=1):
        """
        Evaluate the derivative at parameter.

        :param u:
        :param d:

        :return:
        """
        v = Vector(self.DN(u, d).XYZ())
        return v

    def reverse(self):
        """
        Reverse the direction of the line.

        :return:
        """
        self.Reverse()

    def reversed_u(self, u):
        """
        Calculate the parameter on the reversed line.

        :param u:

        :return:
        """
        return self.ReversedParameter(u)

    def arc_length(self, u1, u2):
        """
        Calculate the curve length between the parameter.

        :param u1:
        :param u2:

        :return:
        """
        return curve_length(self, u1, u2)


class Line(Geom_Line, Curve):
    """
    Line.
    """

    def __init__(self, *args):
        super(Line, self).__init__(*args)
        Curve.__init__(self)


class NurbsCurve(Geom_BSplineCurve, Curve):
    """
    NURBS curve.
    """

    def __init__(self, *args):
        super(NurbsCurve, self).__init__(*args)
        Curve.__init__(self)

    @property
    def p(self):
        return self.Degree()

    @property
    def n(self):
        return self.NbPoles() - 1

    @property
    def knots(self):
        tcol_array = TColStd_Array1OfReal(1, self.NbKnots())
        self.Knots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def mult(self):
        tcol_array = TColStd_Array1OfInteger(1, self.NbKnots())
        self.Multiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbPoles() +
                                             self.Degree() + 1)
        self.KnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def cp(self):
        tcol_array = TColgp_Array1OfPnt(1, self.NbPoles())
        self.Poles(tcol_array)
        return to_np_from_tcolgp_array1_pnt(tcol_array)

    @property
    def w(self):
        tcol_array = TColStd_Array1OfReal(1, self.NbPoles())
        self.Weights(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def cpw(self):
        return homogenize_array1d(self.cp, self.w)

    @property
    def length(self):
        return self.arc_length(self.u1, self.u2)

    def set_domain(self, u1=0., u2=1.):
        """
        Reparameterize the knot vector between *u1* and *u2*.

        :param float u1: Lower domain.
        :param float u2: Upper domain.

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
        Convert parameter(s) from local domain 0 <= u <= 1 to global domain
        a <= u <= b.

        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        return local_to_global_param(self.u1, self.u2, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0 <= u <= 1.

        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        return global_to_local_param(self.u1, self.u2, *args)

    def segment(self, u1, u2):
        """
        Segment the curve between parameters.

        :param u1:
        :param u2:

        :return:
        """
        if u1 > u2:
            return False
        self.Segment(u1, u2)
        return True


class NurbsCurve2D(Geom2d_BSplineCurve, Geometry):
    """
    2-D NURBS curve.
    """

    def __init__(self, *args):
        super(NurbsCurve2D, self).__init__(*args)
        Geometry.__init__(self)

    @property
    def handle(self):
        return self.GetHandle()

    @property
    def adaptor(self):
        return GeomAdaptor_Curve(self.GetHandle())

    @property
    def u1(self):
        return self.FirstParameter()

    @property
    def u2(self):
        return self.LastParameter()

    @property
    def p(self):
        return self.Degree()

    @property
    def n(self):
        return self.NbPoles() - 1

    @property
    def knots(self):
        tcol_array = TColStd_Array1OfReal(1, self.NbKnots())
        self.Knots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def mult(self):
        tcol_array = TColStd_Array1OfInteger(1, self.NbKnots())
        self.Multiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbPoles() +
                                             self.Degree() + 1)
        self.KnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def is_closed(self):
        return self.IsClosed()

    @property
    def is_periodic(self):
        return self.IsPeriodic()

    def set_domain(self, u1=0., u2=1.):
        """
        Reparameterize the knot vector between *u1* and *u2*.

        :param float u1: Lower domain.
        :param float u2: Upper domain.

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
        Convert parameter(s) from local domain 0 <= u <= 1 to global domain
        a <= u <= b.

        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        return local_to_global_param(self.u1, self.u2, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0 <= u <= 1.

        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        return global_to_local_param(self.u1, self.u2, *args)

    def eval(self, u):
        """
        Evaluate curve at parameter.

        :param float u: Curve parameter.

        :return: Point on curve.
        :rtype: :class:`.Point2D` or ndarray
        """
        p = Point2D()
        self.D0(u, p)
        return p

    def reverse(self):
        """
        Reverse the direction of the curve.

        :return:
        """
        self.Reverse()

    def reversed_u(self, u):
        """
        Calculate the parameter on the reversed curve.

        :param u:

        :return:
        """
        return self.ReversedParameter(u)

    def segment(self, u1, u2):
        """
        Segment the curve between parameters.

        :param u1:
        :param u2:

        :return:
        """
        if u1 > u2:
            return False
        self.Segment(u1, u2)
        return True


class Surface(Geom_Surface, Geometry):
    """
    Base class for surfaces.
    """

    def __init__(self):
        Geometry.__init__(self)

    @property
    def handle(self):
        return self.GetHandle()

    @property
    def adaptor(self):
        return GeomAdaptor_Surface(self.GetHandle())

    def eval(self, u=0., v=0.):
        """
        Evaluate plane at parameters.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.

        :return: Point on surface.
        :rtype: :class:`.Point` or ndarray
        """
        p = Point()
        self.D0(u, v, p)
        return p

    def deriv(self, u, v, nu, nv):
        """
        Calculate the plane derivative.

        :param u:
        :param v:
        :param nu:
        :param nv:

        :return:
        """
        v = Vector(self.DN(u, v, nu, nv).XYZ())
        return v

    def norm(self, u, v):
        """
        Calculate the plane normal.

        :param u:
        :param v:

        :return:
        """
        du = self.deriv(u, v, 1, 0)
        dv = self.deriv(u, v, 0, 1)
        vn = Vector(du.Crossed(dv).XYZ())
        return vn

    def __call__(self, u, v, nu=0, nv=0):
        """
        Evaluate surface or its derivative

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.
        :param int nu: derivative order in u-direction
        :param int nv: derivative order in v-direction

        :return: Point on surface or derivative
        :rtype: :class:`.Point` or Vector
        """
        if (nu + nv) < 1:
            return self.eval(u, v)
        else:
            return self.deriv(u, v, nu, nv)

    def eval_params(self, uprms, vprms):
        """
        Evaluate the surface at multiple parameters.

        :param uprms:
        :param vprms:

        :return:
        """
        return [self.eval(u, v) for u, v in zip(uprms, vprms)]


class Plane(Geom_Plane, Surface):
    """
    Plane.
    """

    def __init__(self, *args):
        super(Plane, self).__init__(*args)
        Surface.__init__(self)


class NurbsSurface(Geom_BSplineSurface, Surface):
    """
    NURBS surface.
    """

    def __init__(self, *args):
        super(NurbsSurface, self).__init__(*args)
        Surface.__init__(self)

    @property
    def u1(self):
        return self.Bounds()[0]

    @property
    def u2(self):
        return self.Bounds()[1]

    @property
    def v1(self):
        return self.Bounds()[2]

    @property
    def v2(self):
        return self.Bounds()[3]

    @property
    def p(self):
        return self.UDegree()

    @property
    def q(self):
        return self.VDegree()

    @property
    def n(self):
        return self.NbUPoles() - 1

    @property
    def m(self):
        return self.NbVPoles() - 1

    @property
    def uknots(self):
        tcol_array = TColStd_Array1OfReal(1, self.NbUKnots())
        self.UKnots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def umult(self):
        tcol_array = TColStd_Array1OfInteger(1, self.NbUKnots())
        self.UMultiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def uk(self):
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbUPoles() +
                                             self.UDegree() + 1)
        self.UKnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def vknots(self):
        tcol_array = TColStd_Array1OfReal(1, self.NbVKnots())
        self.VKnots(tcol_array)
        return to_np_from_tcolstd_array1_real(tcol_array)

    @property
    def vmult(self):
        tcol_array = TColStd_Array1OfInteger(1, self.NbVKnots())
        self.VMultiplicities(tcol_array)
        return to_np_from_tcolstd_array1_integer(tcol_array)

    @property
    def vk(self):
        tcol_knot_seq = TColStd_Array1OfReal(1, self.NbVPoles() +
                                             self.VDegree() + 1)
        self.VKnotSequence(tcol_knot_seq)
        return to_np_from_tcolstd_array1_real(tcol_knot_seq)

    @property
    def cp(self):
        tcol_array = TColgp_Array2OfPnt(1, self.NbUPoles(), 1, self.NbVPoles())
        self.Poles(tcol_array)
        return to_np_from_tcolgp_array2_pnt(tcol_array)

    @property
    def w(self):
        tcol_array = TColStd_Array2OfReal(1, self.NbUPoles(),
                                          1, self.NbVPoles())
        self.Weights(tcol_array)
        return to_np_from_tcolstd_array2_real(tcol_array)

    @property
    def cpw(self):
        return homogenize_array2d(self.cp, self.w)

    def set_udomain(self, u1=0., u2=1.):
        """
        Reparameterize the knot vector between *u1* and *u2* in the
        u-direction.

        :param float u1: Lower domain.
        :param float u2: Upper domain.

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
        Reparameterize the knot vector between *v1* and *v2* in the
        v-direction.

        :param float v1: Lower domain.
        :param float v2: Upper domain.

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
        Convert parameter(s) from local domain 0 <= u,v <= 1 to global domain
        a <= u,v <= b.

        :param str d: Direction ('u' or 'v').
        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        if d.lower() in ['u']:
            return local_to_global_param(self.u1, self.u2, *args)
        else:
            return local_to_global_param(self.v1, self.v2, *args)

    def global_to_local_param(self, d, *args):
        """
        Convert surface parameter(s) from global domain a <= u,v <= b to local
        domain 0 <= u,v <= 1.

        :param str d: Direction ('u' or 'v').
        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        if d.lower() in ['u']:
            return global_to_local_param(self.u1, self.u2, *args)
        else:
            return global_to_local_param(self.v1, self.v2, *args)

    def segment(self, u1, u2, v1, v2):
        """
        Segment the surface between the parameters.

        :param u1:
        :param u2:
        :param v1:
        :param v2:

        :return:
        """
        if u1 > u2 or v1 > v2:
            return False
        self.CheckAndSegment(u1, u2, v1, v2)
        return True


# TODO Move to topology.
class BBox(Bnd_Box):
    """
    Bounding box in 3-D.
    """

    def __init__(self):
        super(BBox, self).__init__()

    @property
    def is_void(self):
        return self.IsVoid()

    @property
    def pmin(self):
        if self.is_void:
            return None
        return Point(self.CornerMin().XYZ())

    @property
    def pmax(self):
        if self.is_void:
            return None
        return Point(self.CornerMax().XYZ())

    @property
    def xmin(self):
        if self.is_void:
            return None
        return self.CornerMin().X()

    @property
    def xmax(self):
        if self.is_void:
            return None
        return self.CornerMax().X()

    @property
    def ymin(self):
        if self.is_void:
            return None
        return self.CornerMin().Y()

    @property
    def ymax(self):
        if self.is_void:
            return None
        return self.CornerMax().Y()

    @property
    def zmin(self):
        if self.is_void:
            return None
        return self.CornerMin().Z()

    @property
    def zmax(self):
        if self.is_void:
            return None
        return self.CornerMax().Z()

    def add_box(self, box):
        """
        Add the other box to the box.

        :param box:

        :return:
        """
        if not isinstance(box, Bnd_Box):
            return False
        self.Add(box)
        return True

    def add_shape(self, shape):
        """
        Add shape to the boudning box.

        :param shape:

        :return:
        """
        brepbndlib_Add(shape, self, True)
        return True

    def distance(self, box):
        """
        Calculate distance to other box.

        :param box:

        :return:
        """
        if self.is_void or box.IsVoid():
            return None
        return self.Distance(box)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
