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
__all__ = ["Property", "Shell", "PropertyData"]


class Property(object):
    """
    Base class for properties.
    """
    _all = {}
    _indx = 1

    def __init__(self, name):
        self._name = name
        Property._all[name] = self
        self._pid = Property._indx
        self._export = True
        Property._indx += 1

    def __int__(self):
        return self._pid

    @property
    def name(self):
        return self._name

    @property
    def pid(self):
        return self._pid

    @property
    def export(self):
        return self._export

    @classmethod
    def get_property(cls, prop):
        """
        Get a property.

        :param prop:

        :return:
        """
        if isinstance(prop, Property):
            return prop
        try:
            return Property._all[prop]
        except KeyError:
            return None

    def set_export(self, export=True):
        """
        Set option to export property.

        :param bool export: Option to export property.

        :return: *True* if set to export, *False* if not.
        :rtype: bool
        """
        if export:
            self._export = True
            return True
        self._export = False
        return False


class Shell(Property):
    """
    Shell element property.
    """

    def __init__(self, name, t, mid1):
        super(Shell, self).__init__(name)
        self._t = t
        self._mid1 = mid1

    @property
    def t(self):
        return self._t

    @property
    def mid1(self):
        return int(self._mid1)


class PropertyData(object):
    """
    Property data manager.
    """

    @staticmethod
    def get_property(prop):
        """
        Get a property.

        :param prop:

        :return:
        """
        return Property.get_property(prop)

    @staticmethod
    def set_export(prop, export=True):
        """
        Set option to export material.

        :param prop: Property name or instance.
        :param bool export: Option to export property.

        :return: *True* if set to export, *False* if not.
        :rtype: bool
        """
        prop = PropertyData.get_property(prop)
        if not prop:
            return False
        return prop.set_export(export)

    @staticmethod
    def create_shell(name, t, mid1):
        """
        Create shell element property.

        :param name: Material name.
        :param float t: Default membrane thickness.
        :param mid1: Material for the membrane.

        :return: New shell element property.
        :rtype: :class:`.Shell`
        """
        return Shell(name, t, mid1)
