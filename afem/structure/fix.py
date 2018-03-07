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

from afem.structure.group import Group, GroupAPI
from afem.topology.fix import FixShape

__all__ = ["FixGroup"]


class FixGroup(object):
    """
    Attempt to fix the shapes of each part in an group using
    :class:`.FixShape`. Subgroups are included by default.

    :param group: The group. If ``None`` then the active group is used.
    :type group: str or afem.structure.group.Group or None
    :param float precision: Basic precision value.
    :param float min_tol: Minimum allowed tolerance.
    :param float max_tol: Maximum allowed tolerance.

    :raise TypeError: If an :class:`.Group` instance is not found.
    """

    def __init__(self, group=None, precision=None, min_tol=None, max_tol=None):
        group = GroupAPI.get_group(group)
        if not isinstance(group, Group):
            raise TypeError('Could not find group.')

        parts = group.get_parts()
        compound = group.as_compound()

        fix = FixShape(compound, precision, min_tol, max_tol)

        for part in parts:
            new_shape = fix.apply(part.shape)
            part.set_shape(new_shape)

    @staticmethod
    def limit_tolerance(group=None, tol=1.0e-7):
        """
        Limit tolerances for the group shapes.

        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None
        :param float tol: Target tolerance.

        :return: *True* if at least one tolerance of a sub-shape was modified.
        :rtype: bool

        :raise TypeError: If an :class:`.Group` instance is not found.
        """
        group = GroupAPI.get_group(group)
        if not isinstance(group, Group):
            raise TypeError('Could not find group.')

        shape = group.as_compound()
        return FixShape.limit_tolerance(shape, tol)

    @staticmethod
    def set_tolerance(group=None, tol=1.0e7):
        """
        Enforce tolerance on the given group.

        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None
        :param float tol: The tolerance.

        :return: None.

        :raise TypeError: If an :class:`.Group` instance is not found.
        """
        group = GroupAPI.get_group(group)
        if not isinstance(group, Group):
            raise TypeError('Could not find group.')

        shape = group.as_compound()
        return FixShape.set_tolerance(shape, tol)
