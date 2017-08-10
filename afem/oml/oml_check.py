from afem.oml.oml_entities import Body, Fuselage, Wing

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

    @staticmethod
    def is_wing(entity):
        """
        Check to see if the entity is a wing.

        :param entity: Entity to check.

        :return: *True* if a Wing, *False* if not.
        :rtype: bool
        """
        return isinstance(entity, Wing)

    @staticmethod
    def is_fuselage(entity):
        """
        Check to see if the entity is a fuselage.

        :param entity: Entity to check.

        :return: *True* if a Fuselage, *False* if not.
        :rtype: bool
        """
        return isinstance(entity, Fuselage)
