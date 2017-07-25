from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.BRepGProp import brepgprop
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.GProp import GProp_GProps

from ...config import Settings

__all__ = []


def curve_length(curve, u1, u2):
    """
    Calculate the length of a curve.
    """
    if u1 > u2:
        u1, u2 = u2, u1
    return GCPnts_AbscissaPoint.Length(curve.adaptor, u1, u2, Settings.gtol)


def surface_area(surface):
    """
    Calculate surface area.
    """
    f = BRepBuilderAPI_MakeFace(surface.handle, 0.).Face()
    sprops = GProp_GProps()
    brepgprop.SurfaceProperties(f, sprops, Settings.gtol)
    return sprops.Mass()
