from .checker import CheckPart
from .creator import CreatePart
from .methods.cut_parts import cut_part
from .methods.fuse_parts import fuse_surface_parts, fuse_wing_parts
from .methods.sew_parts import sew_surface_parts
from .methods.split_parts import split_parts


class PartTools(object):
    """
    Part tools.
    """
    create = CreatePart()
    check = CheckPart()

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
        return True in [cut_part(part, cutter) for part in parts if
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

    @staticmethod
    def split_parts(parts, tools=None, order=False):
        """
        Split parts.
        
        :param parts: 
        :param tools:
        :param bool order:
        
        :return: 
        """
        if order:
            parts = PartTools.order_parts(parts)
        return split_parts(parts, tools)

    @staticmethod
    def order_parts(parts):
        """
        Order the list of parts by id.
        
        :param parts:
         
        :return: 
        """
        if not parts:
            return []

        part_order = [(part.id, part) for part in parts]
        part_order.sort(key=lambda tup: tup[0])
        return [row[1] for row in part_order]
