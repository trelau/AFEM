from .fuselage_part import FuselagePart


class Frame(FuselagePart):
    """
    Frame.
    """

    def __init__(self, name, surface_shape):
        super(Frame, self).__init__(name, surface_shape)
