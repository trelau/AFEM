from .methods.join_parts import fuse_wing_parts
from .wing_part import WingPart


class PartTools(object):
    """
    Part tools.
    """

    @staticmethod
    def fuse_wing_parts(parts, tol=None):
        """
        Automatically fuse wing parts.

        :param parts:
        :param float tol:

        :return:
        """
        _wing_parts = []
        for part in parts:
            if isinstance(part, WingPart):
                _wing_parts.append(part)
        return fuse_wing_parts(_wing_parts, tol)
