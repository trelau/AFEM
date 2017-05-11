from OCC.Geom import Geom_Curve, Geom_Surface

from .points import Point


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


Geom_Curve.u1 = property(_u1)
Geom_Curve.u2 = property(_u2)

Geom_Curve.eval = _ceval
Geom_Surface.eval = _seval
