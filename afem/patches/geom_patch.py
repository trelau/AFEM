from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Geom import Geom_Geometry
from OCC.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Standard import Standard_Transient, Handle_Standard_Transient
from OCC.gce import gce_MakeMirror

__all__ = []


def _get_handle(self):
    return self.GetHandle()


def _get_handle_handle(self):
    return self


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


def _set_mirror(self, pln):
    self.mirror = pln


def _get_mirrored(self):
    if not self.mirror:
        return None
    trsf = gce_MakeMirror(self.mirror.Pln()).Value()
    builder = BRepBuilderAPI_Transform(self, trsf, True)
    if not builder.IsDone():
        return None
    return builder.Shape()


Standard_Transient.handle = property(_get_handle)
Handle_Standard_Transient.handle = property(_get_handle_handle)

Geom_Geometry.color = None
Geom_Geometry.set_color = _set_color
Geom_Geometry.transparency = 0.
Geom_Geometry.set_transparency = _set_transparency
Geom_Geometry.mirror = None
Geom_Geometry.set_mirror = _set_mirror
Geom_Geometry.get_mirrored = _get_mirrored
