from .curve_part import CurvePart
from .surface_part import SurfacePart


class Stiffener1D(CurvePart):
    """
    1-D stiffener for surface parts.
    """

    def __init__(self, label, shape, cref=None):
        super(Stiffener1D, self).__init__(label, shape, cref, False)


class Stiffener2D(SurfacePart):
    """
    2-D stiffener for surface parts.
    """

    def __init__(self, label, shape, cref=None, sref=None):
        super(Stiffener2D, self).__init__(label, shape, cref, sref, False)
