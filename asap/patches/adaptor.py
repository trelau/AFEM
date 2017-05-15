from OCC.Adaptor3d import Adaptor3d_Curve, Adaptor3d_Surface
from OCC.GCPnts import GCPnts_AbscissaPoint

from ..geometry.points import Point


def _u1(self):
    """
    First parameter of curve.
    """
    return self.FirstParameter()


def _u2(self):
    """
    Last parameter of curve.
    """
    return self.LastParameter()


def _ceval(self, u):
    """
    Evaluate a point on the curve.
    """
    p = Point()
    self.D0(u, p)
    return p


def _seval(self, u, v):
    """
    Evaluate a point on the surface.
    """
    p = Point()
    self.D0(u, v, p)
    return p


def _abscissa_point(self, dx, u0=None):
    """
    Evaluate point on curve.
    """
    if u0 is None:
        u0 = self.FirstParameter()
    elif u0 < self.FirstParameter():
        u0 = self.FirstParameter()
    elif u0 > self.LastParameter():
        u0 = self.LastParameter()
    tool = GCPnts_AbscissaPoint(self, dx, u0)
    if not tool.IsDone():
        return None
    u = tool.Parameter()
    p = Point()
    self.D0(u, p)
    return p


Adaptor3d_Curve.u1 = property(_u1)
Adaptor3d_Curve.u2 = property(_u2)
Adaptor3d_Curve.eval = _ceval
Adaptor3d_Curve.eval_dx = _abscissa_point

Adaptor3d_Surface.eval = _seval
