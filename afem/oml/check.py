from afem.oml.entities import Body

__all__ = ["CheckOML"]


class CheckOML(object):
    """
    Check OML.
    """

    @staticmethod
    def is_body(entity):
        """
        Check to see if the entity is a Body.

        :param entity: Entity to check.

        :return: *True* if Body, *False* if not.
        :rtype: bool
        """
        return isinstance(entity, Body)
