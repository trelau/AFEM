from .fuselage_part import FuselagePart


class Floor(FuselagePart):
    """
    Floor.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Floor, self).__init__(label, shape, cref, sref, assy)
