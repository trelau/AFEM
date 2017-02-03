from .wing_part import WingPart


class Spar(WingPart):
    """
    Wing spar.
    """

    def __init__(self, name, wing, rshape):
        super(Spar, self).__init__(name, wing, rshape)
