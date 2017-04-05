from .wing_part import WingPart


class Rib(WingPart):
    """
    Wing rib.
    """

    def __init__(self, name, wing, surface_shape):
        super(Rib, self).__init__(name, wing, surface_shape)
