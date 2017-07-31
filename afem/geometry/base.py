from ..graphics.viewer import ViewableItem

__all__ = ["Geometry"]


class Geometry(ViewableItem):
    """
    Base class for geometry.
    """

    def __init__(self):
        super(Geometry, self).__init__()
