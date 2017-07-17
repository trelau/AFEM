from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, \
    BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.GCPnts import GCPnts_AbscissaPoint
from OCC.Geom import Geom_Curve, Geom_Geometry, Geom_Surface
from OCC.GeomAdaptor import GeomAdaptor_Curve, GeomAdaptor_Surface
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Standard import Standard_Transient
from OCC.gp import gp_Pnt
from numpy import array, float64

from ..config import Settings
from ..geometry.points import Point
from ..geometry.projector import ProjectGeom


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


def _get_handle(self):
    return self.GetHandle()


def _curve_adaptor(self):
    return GeomAdaptor_Curve(self.GetHandle())


def _surface_adaptor(self):
    return GeomAdaptor_Surface(self.GetHandle())


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


def _edge(self):
    """
    Create edge from curve.
    """
    return BRepBuilderAPI_MakeEdge(self.GetHandle()).Edge()


def _to_edge(self, u1, u2):
    """
    Make edge from curve between parameters or points.
    """
    return BRepBuilderAPI_MakeEdge(self.GetHandle(), u1, u2).Edge()


def _wire(self):
    """
    Make wire from curve.
    """
    e = BRepBuilderAPI_MakeEdge(self.GetHandle()).Edge()
    return BRepBuilderAPI_MakeWire(e).Wire()


def _to_wire(self, u1, u2):
    """
    Make wire from curve between parameters or points.
    """
    e = BRepBuilderAPI_MakeEdge(self.GetHandle(), u1, u2).Edge()
    return BRepBuilderAPI_MakeWire(e).Wire()


def _face(self):
    """
    Make face from surface. 
    """
    return BRepBuilderAPI_MakeFace(self.GetHandle()).Face()


def _to_face(self, u1, u2, v1, v2):
    """
    Make face from surface trimmed by parameters.
    """
    return BRepBuilderAPI_MakeFace(self.GetHandle(), u1, u2, v1, v2).Face()


def _length(self):
    """
    Curve arc length.
    """
    return GCPnts_AbscissaPoint.Length(GeomAdaptor_Curve(self.GetHandle()),
                                       Settings.gtol)


def _xyz(self):
    """
    Coordinates of gp_Pnt.
    """
    return array([self.X(), self.Y(), self.Z()], dtype=float64)


def _array(self, dtype=float64, copy=True, order=None, subok=False,
           ndmin=0):
    return array(_xyz(self), dtype=dtype, copy=copy, order=order, subok=subok,
                 ndmin=ndmin)


def _string(self):
    """
    Print gp_pnt coordinates.
    """
    return 'gp_Pnt = ({0}, {1}, {2})'.format(*self.xyz)


def _invert_pnt_on_geom(self, pnt):
    """
    Find parameter of point on geometry.
    """
    return ProjectGeom.invert(pnt, self)


def _project_pnt_to_geom(self, pnt):
    """
    Project point to geometry.
    """
    return ProjectGeom.point_to_geom(pnt, self, True)


def _abscissa_point(self, dx, u0=None, is_local=False):
    """
    Evaluate point on curve.
    """
    adp_crv = _curve_adaptor(self)
    if u0 is None:
        u0 = adp_crv.FirstParameter()
    elif u0 < adp_crv.FirstParameter():
        u0 = adp_crv.FirstParameter()
    elif u0 > adp_crv.LastParameter():
        u0 = adp_crv.LastParameter()
    # Multiply dx by length if is_local=True
    if is_local:
        dx = dx * _length(self)
    tool = GCPnts_AbscissaPoint(adp_crv, dx, u0)
    if not tool.IsDone():
        return None
    u = tool.Parameter()
    p = Point()
    adp_crv.D0(u, p)
    return p


def _set_color(self, r, g, b):
    """
    Set color of shape.
    """
    if r > 1.:
        r /= 255.
    if g > 1.:
        g /= 255.
    if b > 1.:
        b /= 255.
    self.color = Quantity_Color(r, g, b, Quantity_TOC_RGB)


def _set_transparency(self, transparency):
    """
    Set transparency of shape.
    """
    self.transparency = transparency


Standard_Transient.handle = property(_get_handle)

Geom_Geometry.color = None
Geom_Geometry.set_color = _set_color
Geom_Geometry.transparency = 0.
Geom_Geometry.set_transparency = _set_transparency

Geom_Curve.u1 = property(_u1)
Geom_Curve.u2 = property(_u2)
Geom_Curve.length = property(_length)
Geom_Curve.adaptor = property(_curve_adaptor)
Geom_Curve.edge = property(_edge)
Geom_Curve.to_edge = _to_edge
Geom_Curve.wire = property(_wire)
Geom_Curve.to_wire = _to_wire
Geom_Curve.eval = _ceval
Geom_Curve.invert = _invert_pnt_on_geom
Geom_Curve.project = _project_pnt_to_geom
Geom_Curve.eval_dx = _abscissa_point

Geom_Surface.face = property(_face)
Geom_Surface.to_face = _to_face
Geom_Surface.adaptor = property(_surface_adaptor)
Geom_Surface.eval = _seval
Geom_Surface.invert = _invert_pnt_on_geom
Geom_Surface.project = _project_pnt_to_geom

gp_Pnt.xyz = property(_xyz)
gp_Pnt.__array__ = _array
gp_Pnt.__str__ = _string
