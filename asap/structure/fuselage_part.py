from .surface_part import SurfacePart


class FuselagePart(SurfacePart):
    """
    Base class for fuselage parts.
    """

    def __init__(self, name, fuselage, rshape):
        super(FuselagePart, self).__init__(name, rshape)
        self._fuselage = fuselage
