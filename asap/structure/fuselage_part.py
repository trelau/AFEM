from .surface_part import SurfacePart


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(FuselagePart, self).__init__(label, shape, cref, sref)
