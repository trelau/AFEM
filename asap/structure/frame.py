from .fuselage_part import FuselagePart


class Frame(FuselagePart):
    """
    Frame.
    """

    def __init__(self, name, fuselage, surface_shape):
        super(Frame, self).__init__(name, fuselage, surface_shape)
