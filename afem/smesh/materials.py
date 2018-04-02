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
__all__ = ["Material", "Isotropic", "MaterialData"]


class Material(object):
    """
    Base class for materials.
    """
    _all = {}
    _indx = 1

    def __init__(self, name):
        self._name = name
        Material._all[name] = self
        self._mid = Material._indx
        self._export = True
        Material._indx += 1

    def __int__(self):
        return self._mid

    @property
    def name(self):
        return self._name

    @property
    def mid(self):
        return self._mid

    @property
    def export(self):
        return self._export

    @classmethod
    def get_material(cls, material):
        """
        Get a material.

        :param material:

        :return:
        """
        if isinstance(material, Material):
            return material
        try:
            return Material._all[material]
        except KeyError:
            return None

    def set_export(self, export=True):
        """
        Set option to export material.

        :param bool export: Option to export material.

        :return: *True* if set to export, *False* if not.
        :rtype: bool
        """
        if export:
            self._export = True
            return True
        self._export = False
        return False


class Isotropic(Material):
    """
    Linear isotropic material.
    """

    def __init__(self, name, E, G, nu, rho):
        super(Isotropic, self).__init__(name)
        self._E = E
        self._G = G
        self._nu = nu
        self._rho = rho

    @property
    def E(self):
        return self._E

    @property
    def G(self):
        return self._G

    @property
    def nu(self):
        return self._nu

    @property
    def rho(self):
        return self._rho


class MaterialData(object):
    """
    Material data manager.
    """

    @staticmethod
    def get_material(material):
        """
        Get a material.

        :param material:

        :return:
        """
        return Material.get_material(material)

    @staticmethod
    def set_export(material, export=True):
        """
        Set option to export material.

        :param material: Material name or instance.
        :param bool export: Option to export material.

        :return: *True* if set to export, *False* if not.
        :rtype: bool
        """
        mat = MaterialData.get_material(material)
        if not mat:
            return False
        return mat.set_export(export)

    @staticmethod
    def create_isotropic(name, E, G, nu, rho):
        """
        Create linear isotropic material.

        :param name: Material name.
        :param E: Elastic modulus.
        :param G: Shear modulus.
        :param nu: Poisson's ratio.
        :param rho: Density.

        :return: New linear isotropic material.
        :rtype: :class:`.Isotropic`
        """
        return Isotropic(name, E, G, nu, rho)
