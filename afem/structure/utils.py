# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
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
