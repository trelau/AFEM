from .wing_part import WingPart


class Rib(WingPart):
    """
    Wing rib.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Rib, self).__init__(label, shape, cref, sref, assy)
