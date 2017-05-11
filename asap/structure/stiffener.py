from .curve_part import CurvePart


class Stiffener(CurvePart):
    """
    Stiffener for surface parts.
    """

    def __init__(self, label, shape, cref=None):
        super(Stiffener, self).__init__(label, shape, cref, False)
