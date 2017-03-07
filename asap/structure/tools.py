from .checker import CheckPart
from .methods.cut_parts import cut_surface_part
from .methods.fuse_parts import fuse_surface_parts, fuse_wing_parts
from .methods.sew_parts import sew_surface_parts


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
            if CheckPart.is_wing_part(part):
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
        return True in [cut_surface_part(part, cutter) for part in parts if
                        CheckPart.is_surface_part(part)]

    @staticmethod
    def sew_parts(parts):
        """
        Sew surface parts.

        :param parts:

        :return:
        """
        _parts = [part for part in parts if CheckPart.is_part(part)]
        return sew_surface_parts(_parts)

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
