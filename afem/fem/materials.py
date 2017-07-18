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
