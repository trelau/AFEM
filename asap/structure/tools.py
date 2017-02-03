from .methods.join_parts import join_wing_parts
from .wing_part import WingPart


class PartTools(object):
    """
    Part tools.
    """

    @staticmethod
    def join_wing_parts(parts, tol=None):
        """
        Automatically join wing parts.

        :param parts:
        :param float tol:

        :return:
        """
        _wing_parts = []
        for part in parts:
            if isinstance(part, WingPart):
                _wing_parts.append(part)
        return join_wing_parts(_wing_parts, tol)
