from .part import Part
from ..topology import ShapeTools


class CurvePart(Part):
    """
    Base class for curve-based parts.
    """

    def __init__(self, label, shape, cref=None, add_to_assy=True):
        super(CurvePart, self).__init__(label, shape, cref, None, add_to_assy)

    @property
    def edges(self):
        return ShapeTools.get_edges(self)

    @property
    def nedges(self):
        return len(self.edges)

    @property
    def reshapes(self):
        return self.edges
