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
from afem.base.entities import NamedItem

__all__ = ["Property", "Shell"]


class Property(NamedItem):
    """
    Base class for properties.
    """

    def __init__(self, name):
        super(Property, self).__init__(name)


class Shell(Property):
    """
    Shell element property.
    """

    def __init__(self, name, t):
        super(Shell, self).__init__(name)
        self._t = t

    @property
    def t(self):
        return self._t
