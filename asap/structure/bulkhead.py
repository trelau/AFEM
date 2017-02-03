from .fuselage_part import FuselagePart


class Bulkhead(FuselagePart):
    """
    Bulkhead.
    """

    def __init__(self, name, fuselage, rshape):
        super(Bulkhead, self).__init__(name, fuselage, rshape)
