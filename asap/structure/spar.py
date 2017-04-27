from .wing_part import WingPart


class Spar(WingPart):
    """
    Wing spar.
    """

    def __init__(self, label, surface_shape):
        super(Spar, self).__init__(label, surface_shape)
