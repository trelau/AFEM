from .wing_part import WingPart


class Rib(WingPart):
    """
    Wing rib.
    """

    def __init__(self, label, shape, cref=None, sref=None, assy=None):
        super(Rib, self).__init__(label, shape, cref, sref, assy)


class RibByPoints(Rib):
    """
    Rib between two points.
    """

    def __init__(self, label, wing, p1, p2):
        # Insert method to find shape here
        shape = None
        super(RibByPoints, self).__init__(label, shape)
