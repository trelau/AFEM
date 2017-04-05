from .fuselage_part import FuselagePart


class Floor(FuselagePart):
    """
    Floor.
    """

    def __init__(self, name, fuselage, surface_shape):
        super(Floor, self).__init__(name, fuselage, surface_shape)
