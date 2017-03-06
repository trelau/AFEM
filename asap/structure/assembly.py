class Assembly(object):
    """
    Structural assembly.
    """

    def __init__(self, name, parent):
        self._name = name
        self._parent = parent
        self._children = set()
        self._parts = set()

    @property
    def name(self):
        return self._name

    @property
    def parent(self):
        return self._parent

    @property
    def parts(self):
        return list(self._parts)

    @property
    def elements(self):
        elms = {}
        for part in self.parts:
            for e in part.elements:
                elms[e.eid] = e
        return elms.values()

    @property
    def nodes(self):
        nodes = {}
        for part in self.parts:
            for n in part.nodes:
                nodes[n.nid] = n
        return nodes.values()

    def activate(self):
        """
        Activate this assembly.

        :return:
        """
        AssemblyData._active = self

    def add_parts(self, *parts):
        """
        Add parts to the assembly.

        :param parts: Parts to add.
        """
        part_set = set(parts)
        self._parts.update(part_set)

    def get_part(self, name):
        """
        Get a part in the assembly by name.

        :param name: Name of part.

        :return: Part or *None* if part is not found.
        """
        for part in self.parts:
            if part.name == name:
                return part
        return None

    def mesh(self, maxh=1., quad_dominated=True):
        """
        Mesh the assembly.

        :param float maxh:
        :param bool quad_dominated: Option to generate a quad-dominated mesh.

        :return: Mesh status dictionary where the part names are the key
            and the status is the value.
        :rtype: dict
        """
        # Mesh parts in assembly.
        results = {}
        for part in self.parts:
            status = part.mesh(maxh, quad_dominated=quad_dominated)
            results[part.name] = status

        return results


class AssemblyData(object):
    """
    Assembly data.
    """
    _master = Assembly('Model', None)
    _all = {'Model': _master}
    _active = _master

    @classmethod
    def get_active(cls):
        """
        Get the active assembly.

        :return: Active assembly.
        :rtype: :class:`.Assembly`
        """
        assy = cls._active
        if not isinstance(assy, Assembly):
            return cls._master
        return assy

    @classmethod
    def make_active(cls, assy):
        """
        Activate the assembly.

        :param assy:
        :return:
        """
        assy = cls.get_assy(assy)
        assy.activate()

    @classmethod
    def get_assy(cls, assy=None):
        """
        Get assembly.

        :param assy: Assembly name to get. If an assembly instance is
            provided, then it is returned.
        :return: Assembly or active assembly if name does not exist.
        :rtype: :class:`.Assembly`
        """
        if isinstance(assy, Assembly):
            return assy
        try:
            return cls._all[assy]
        except KeyError:
            return cls.get_active()

    @classmethod
    def create_assy(cls, name, parent=None, active=True, *parts):
        """
        Create an assembly.

        :param name:
        :param parent:
        :param active:
        :param parts:

        :return:
        """
        if name in cls._all:
            return None
        parent = cls.get_assy(parent)
        assy = Assembly(name, parent)
        if active:
            cls._active = assy
        assy.add_parts(*parts)
        cls._all[name] = assy
        return assy

    @classmethod
    def add_parts(cls, assy, *parts):
        """
        Add parts to the assembly.

        :param assy:
        :param parts:

        :return:
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return False
        assy.add_parts(*parts)
        return True

    @classmethod
    def get_parts(cls, assy=None):
        """
        Get parts from assembly.

        :param assy:

        :return:
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return False
        return assy.parts

    @classmethod
    def mesh_assy(cls, assy=None, maxh=1., quad_dominated=True):
        """
        Mesh the assembly.


        :param assy:
        :param float maxh:
        :param bool quad_dominated: Option to generate a quad-dominated mesh.

        :return:
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return {}
        return assy.mesh(maxh, quad_dominated=quad_dominated)
