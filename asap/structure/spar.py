from .wing_part import WingPart


class Spar(WingPart):
    """
    Wing spar.
    """

    def __init__(self, name, wing, surface_shape):
        super(Spar, self).__init__(name, wing, surface_shape)
