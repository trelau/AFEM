from OCC.gp import gp_Ax1, gp_Ax3


class Axis1(gp_Ax1):
    """
    Coordinate system.
    """

    def __init__(self, *args):
        super(Axis1, self).__init__(*args)


class Axis3(gp_Ax3):
    """
    Coordinate system.
    """

    def __init__(self, *args):
        super(Axis3, self).__init__(*args)
