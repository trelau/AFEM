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
