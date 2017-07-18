from .wing_part import WingPart


class Spar(WingPart):
    """
    Wing spar.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Spar, self).__init__(label, shape, cref, sref, assy)
