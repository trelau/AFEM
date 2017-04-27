from .wing_part import WingPart


class Rib(WingPart):
    """
    Wing rib.
    """

    def __init__(self, label, surface_shape):
        super(Rib, self).__init__(label, surface_shape)
