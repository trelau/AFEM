from .body import Body
from .fuselage import Fuselage
from .wing import Wing


class CheckOML(object):
    """
    Check OML components.
    """

    @staticmethod
    def is_body(oml):
        """
        Check to see if the entity is a Body.

        :param oml: Object to check.

        :return: *True* if Body, *False* if not.
        :rtype: bool
        """
        return isinstance(oml, Body)

    @staticmethod
    def is_wing(oml):
        """
        Check to see if the OML component is a wing.

        :param oml: Object to check.

        :return: *True* if a wing, *False* if not.
        :rtype: bool
        """
        return isinstance(oml, Wing)

    @staticmethod
    def is_fuselage(oml):
        """
        Check to see if the OML component is a fuselage.

        :param oml: Object to check.

        :return: *True* if a fuselage, *False* if not.
        :rtype: bool
        """
        return isinstance(oml, Fuselage)
