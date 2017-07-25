from .surface_part import SurfacePart

__all__ = ["Stringer"]


class Stringer(SurfacePart):
    """
    Discretely modeled stringer.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Stringer, self).__init__(label, shape, cref, sref, assy)
