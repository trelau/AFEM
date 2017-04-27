from .fuselage_part import FuselagePart


class Floor(FuselagePart):
    """
    Floor.
    """

    def __init__(self, label, surface_shape):
        super(Floor, self).__init__(label, surface_shape)
