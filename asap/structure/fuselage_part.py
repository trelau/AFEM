from .surface_part import SurfacePart


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """

    def __init__(self, name, fuselage, surface_shape):
        super(FuselagePart, self).__init__(name, surface_shape)
        self._fuselage = fuselage
