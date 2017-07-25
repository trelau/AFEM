from .body import Body

__all__ = ["Fuselage"]


class Fuselage(Body):
    """
    Fuselage body.
    """

    def __init__(self, shape):
        super(Fuselage, self).__init__(shape)
