from .part import Part
from ..topology import ShapeTools

__all__ = ["CurvePart"]


class CurvePart(Part):
    """
    Base class for curve-based parts.
    """

    def __init__(self, label, shape, cref=None, assy=None):
        super(CurvePart, self).__init__(label, shape, cref, None, assy)

    @property
    def reshapes(self):
        return self.edges

    @property
    def length(self):
        return ShapeTools.shape_length(self)
