from .fuselage_part import FuselagePart


class Frame(FuselagePart):
    """
    Frame.
    """

    def __init__(self, name, shape, cref=None, sref=None, assy=None):
        super(Frame, self).__init__(name, shape, cref, sref, assy)
