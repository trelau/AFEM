from OCC.Geom import Geom_BSplineCurve, Geom_Line
from OCC.Geom2d import Geom2d_BSplineCurve
from OCC.GeomAdaptor import GeomAdaptor_Curve
from OCC.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal
from OCC.TColgp import TColgp_Array1OfPnt

from .geom import Geometry
from .methods.calculate import curve_length
from .methods.geom_utils import global_to_local_param, homogenize_array1d, \
    local_to_global_param
from .methods.parameterize import reparameterize_knots
from .points import Point, Point2D
from .vectors import Vector
from ..utils.tcol import to_np_from_tcolgp_array1_pnt, \
    to_np_from_tcolstd_array1_integer, to_np_from_tcolstd_array1_real

__all__ = ["Line", "NurbsCurve", "NurbsCurve2D"]


class Line(Geom_Line, Geometry):
    """
    Line.
    """

    def __init__(self, *args):
        super(Line, self).__init__(*args)
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


class NurbsCurve(Geom_BSplineCurve, Geometry):
    """
    NURBS curve.
    """

    def __init__(self, *args):
        super(NurbsCurve, self).__init__(*args)
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
    def is_closed(self):
        return self.IsClosed()

    @property
    def is_periodic(self):
        return self.IsPeriodic()

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

    def eval(self, u):
        """
        Evaluate curve at parameter.

        :param float u: Curve parameter.

        :return: Point on curve.
        :rtype: :class:`.Point` or ndarray
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
        return Vector(self.DN(u, d).XYZ())

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

    def arc_length(self, u1, u2):
        """
        Calculate the curve length between the parameter.

        :param u1:
        :param u2:

        :return:
        """
        return curve_length(self, u1, u2)

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
