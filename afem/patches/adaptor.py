from OCC.Adaptor3d import Adaptor3d_Curve, Adaptor3d_Surface
from OCC.GCPnts import GCPnts_AbscissaPoint

from ..config import Settings
from ..geometry.points import Point

__all__ = []


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


def _length(self):
    """
    Curve arc length.
    """
    return GCPnts_AbscissaPoint.Length(self, Settings.gtol)


def _abscissa_point(self, dx, u0=None, is_local=False):
    """
    Evaluate point on curve.
    """
    if u0 is None:
        u0 = self.FirstParameter()
    elif u0 < self.FirstParameter():
        u0 = self.FirstParameter()
    elif u0 > self.LastParameter():
        u0 = self.LastParameter()
    # Multiply dx by length if is_local=True
    if is_local:
        dx = dx * _length(self)
    tool = GCPnts_AbscissaPoint(self, dx, u0)
    if not tool.IsDone():
        return None
    u = tool.Parameter()
    p = Point()
    self.D0(u, p)
    return p


Adaptor3d_Curve.u1 = property(_u1)
Adaptor3d_Curve.u2 = property(_u2)
Adaptor3d_Curve.length = property(_length)
Adaptor3d_Curve.eval = _ceval
Adaptor3d_Curve.eval_dx = _abscissa_point

Adaptor3d_Surface.eval = _seval
