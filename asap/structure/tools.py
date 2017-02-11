from .checker import CheckPart
from .methods.fuse_parts import fuse_surface_parts, fuse_wing_parts
from .wing_part import WingPart


class PartTools(object):
    """
    Part tools.
    """

    @staticmethod
    def fuse_parts(parts):
        """
        Fuse surface parts.

        :param parts:

        :return:
        """
        _parts = [part for part in parts if CheckPart.is_part(part)]
        return fuse_surface_parts(_parts)

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
