from .surface_part import SurfacePart

__all__ = ["FuselagePart"]


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(FuselagePart, self).__init__(label, shape, cref, sref, assy)
