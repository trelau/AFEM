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

__all__ = ["order_parts_by_id"]


def order_parts_by_id(parts):
    """
    Order the list of parts by id.

    :param list[afem.structure.entities.Part] parts: The parts.

    :return: List of parts sorted by their ID.
    :rtype: list[afem.structure.entities.Part]
    """
    part_order = [(part.id, part) for part in parts]
    part_order.sort(key=lambda tup: tup[0])
    return [row[1] for row in part_order]
