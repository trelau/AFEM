from .part import Part
from ..topology import ShapeTools


class CurvePart(Part):
    """
    Base class for curve-based parts.
    """

    def __init__(self, label, shape, cref=None, add_to_assy=True):
        super(CurvePart, self).__init__(label, shape, cref, None, add_to_assy)

    @property
    def reshapes(self):
        return self.edges
