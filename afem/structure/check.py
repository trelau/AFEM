#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

from afem.structure.entities import *

__all__ = ["CheckPart"]


class CheckPart(object):
    """
    Check structural components.
    """

    @staticmethod
    def is_part(part):
        """
        Check to see if the part is a structural component.

        :param part: Part to check.

        :return: *True* if part is a structural component, *False* if not.
        :rtype: bool
        """
        return isinstance(part, Part)

    @staticmethod
    def are_parts(*parts):
        """
        Check to see if all entities are parts.

        :param parts: Parts to check.

        :return: *True* if all entities are parts, *False* if not.
        :rtype: bool
        """
        for part in parts:
            if not isinstance(part, Part):
                return False
        return True

    @staticmethod
    def is_surface_part(part):
        """
        Check to see if part is surface-based part.

        :param part:
        :return:
        """
        return isinstance(part, SurfacePart)

    @staticmethod
    def is_wing_part(part):
        """
        Check to see if the part is a wing structural component.

        :param part: Part to check.

        :return: *True* if part is a wing structural component, *False* if
            not.
        :rtype: bool
        """
        return isinstance(part, WingPart)

    @staticmethod
    def is_spar(part):
        """
        Check to see if the part is a spar.

        :param part: Part to check.

        :return: *True* if the part is a spar, *False* if not.
        :rtype: bool
        """
        return isinstance(part, Spar)

    @staticmethod
    def is_rib(part):
        """
        Check to see if the part is a rib.

        :param part: Part to check.

        :return: *True* if the part is a rib, *False* if not.
        :rtype: bool
        """
        return isinstance(part, Rib)
