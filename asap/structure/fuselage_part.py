from .surface_part import SurfacePart


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """

    def __init__(self, label, surface_shape):
        super(FuselagePart, self).__init__(label, surface_shape)
