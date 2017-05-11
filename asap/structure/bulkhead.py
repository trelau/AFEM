from .fuselage_part import FuselagePart


class Bulkhead(FuselagePart):
    """
    Bulkhead.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Bulkhead, self).__init__(label, shape, cref, sref)
