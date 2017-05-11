from .surface_part import SurfacePart


class Stringer(SurfacePart):
    """
    Discretely modeled stringer.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Stringer, self).__init__(label, shape, cref, sref)
