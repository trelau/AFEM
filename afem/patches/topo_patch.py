from OCCT.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCCT.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCCT.ShapeAnalysis import ShapeAnalysis_Edge
from OCCT.TopoDS import TopoDS_Shape
from OCCT.gce import gce_MakeMirror

from afem.topology.check import CheckShape

_analysis_edge = ShapeAnalysis_Edge()

__all__ = []


def _shape(self):
    return self


def _to_shape(self):
    return CheckShape.to_shape(self)


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
    trsf = gce_MakeMirror(self.mirror.handle.Pln()).Value()
    builder = BRepBuilderAPI_Transform(self, trsf, True)
    if not builder.IsDone():
        return None
    return builder.Shape()


TopoDS_Shape.shape = property(_shape)
TopoDS_Shape.downcast = property(_to_shape)

TopoDS_Shape.color = None
TopoDS_Shape.set_color = _set_color
TopoDS_Shape.transparency = 0.
TopoDS_Shape.set_transparency = _set_transparency
TopoDS_Shape.mirror = None
TopoDS_Shape.set_mirror = _set_mirror
TopoDS_Shape.get_mirrored = _get_mirrored
