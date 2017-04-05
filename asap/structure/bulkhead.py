from .fuselage_part import FuselagePart


class Bulkhead(FuselagePart):
    """
    Bulkhead.
    """

    def __init__(self, name, fuselage, surface_shape):
        super(Bulkhead, self).__init__(name, fuselage, surface_shape)
