from OCC.Geom import Geom_BSplineSurface, Geom_Plane
from OCC.GeomAdaptor import GeomAdaptor_Surface
from OCC.TColStd import TColStd_Array1OfInteger, TColStd_Array1OfReal, \
    TColStd_Array2OfReal
from OCC.TColgp import TColgp_Array2OfPnt

from .geom import Geometry
from .methods.geom_utils import global_to_local_param, homogenize_array2d, \
    local_to_global_param
from .methods.parameterize import reparameterize_knots
from .points import Point
from .vectors import Vector
from ..utils.tcol import to_np_from_tcolgp_array2_pnt, \
    to_np_from_tcolstd_array1_integer, to_np_from_tcolstd_array1_real, \
    to_np_from_tcolstd_array2_real


class Plane(Geom_Plane, Geometry):
    """
    Plane.
    """

    def __init__(self, *args):
        super(Plane, self).__init__(*args)
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


class NurbsSurface(Geom_BSplineSurface, Geometry):
    """
    NURBS surface.
    """

    def __init__(self, *args):
        super(NurbsSurface, self).__init__(*args)
        Geometry.__init__(self)

    @property
    def handle(self):
        return self.GetHandle()

    @property
    def adaptor(self):
        return GeomAdaptor_Surface(self.GetHandle())

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

    def D0(self, u, v, p):
        """
        wrap superclass D0 to accept float-like inputs.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.
        :param Point p

        :return:
        :rtype: :class: None
        """
        super(NurbsSurface, self).D0(float(u), float(v), p)
        return None

    def DN(self, u, v, nu, nv):
        """
        wrap superclass DN to accept float-like inputs

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.
        :param int nu: derivative order in u-direction
        :param int nv: derivative order in v-direction

        :return:
        :rtype: :class: gp_Vec
        """
        return super(NurbsSurface, self).DN(float(u), float(v), nu, nv)

    def eval(self, u, v):
        """
        Evaluate surface at parameters.

        :param float u: Parameter in u-direction.
        :param float v: Parameter in v-direction.

        :return: Point on surface.
        :rtype: :class:`.Point` or ndarray
        """
        p = Point()
        self.D0(u, v, p)
        return p

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

    def deriv(self, u, v, nu, nv):
        """
        Calculate the surface derivative.

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
        Calculate the surface normal.

        :param u:
        :param v:

        :return:
        """
        du = self.deriv(u, v, 1, 0)
        dv = self.deriv(u, v, 0, 1)
        vn = Vector(du.Crossed(dv).XYZ())
        return vn

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
