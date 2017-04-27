from .fuselage_part import FuselagePart


class Bulkhead(FuselagePart):
    """
    Bulkhead.
    """

    def __init__(self, label, surface_shape):
        super(Bulkhead, self).__init__(label, surface_shape)
