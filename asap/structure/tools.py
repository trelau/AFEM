from .checker import CheckPart
from .methods.cut_parts import cut_surface_parts
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

    @staticmethod
    def cut_parts(parts, cutter):
        """
        Cut surface parts.

        :param parts:
        :param cutter:

        :return:
        """
        _parts = [part for part in parts if CheckPart.is_part(part)]
        return cut_surface_parts(_parts, cutter)

    @staticmethod
    def discard_faces(parts, shape=None, tol=None):
        """
        Discard faces of parts.

        :param parts:
        :param shape:
        :param tol:

        :return:
        """
        _parts = [part for part in parts if CheckPart.is_surface_part(part)]
        for part in _parts:
            part.discard(shape, tol)
        return True
