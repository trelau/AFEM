from afem.topology.topo_create import CompoundByShapes

__all__ = ["Assembly", "AssemblyData"]


class Assembly(object):
    """
    Structural assembly.
    """

    def __init__(self, label, parent):
        self._label = label
        self._parent = parent
        self._children = set()
        self._parts = set()
        if isinstance(self._parent, Assembly):
            self._parent._children.add(self)
        self._metadata = {}

    @property
    def label(self):
        return self._label

    @property
    def parent(self):
        return self._parent

    @property
    def parts(self):
        return list(self._parts)

    @property
    def elements(self):
        elm_set = set()
        for part in self.parts:
            elm_set.update(part.elements)
        return elm_set

    @property
    def nodes(self):
        node_set = set()
        for part in self.parts:
            node_set.update(part.nodes)
        return node_set

    @property
    def metadata(self):
        return self._metadata

    def add_metadata(self, key, value):
        """
        Add metadata to the assembly.

        :param key:
        :param value:

        :return:
        """
        self._metadata[key] = value
        return True

    def get_metadata(self, key):
        """
        Get metadata.

        :param key:

        :return:
        """
        try:
            return self._metadata[key]
        except KeyError:
            return None

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

    def get_part(self, label):
        """
        Get a part in the assembly by name.

        :param label: Part label.

        :return: Part or *None* if part is not found.
        """
        for part in self.parts:
            if part.label == label:
                return part
        return None

    def get_parts(self, include_subassy=True):
        """
        Get all the parts from the assembly and its sub-assemblies.

        :param include_subassy:

        :return:
        """
        parts = self.parts
        if include_subassy:
            for assy in self._children:
                parts += assy.get_parts(True)
        return parts

    def prepare_shape_to_mesh(self, include_subassy=True):
        """
        Prepare a shape to mesh using the parts in the assembly and its
        sub-assemblies.

        :param include_subassy:

        :return:
        """
        parts = self.get_parts(include_subassy)
        if not parts:
            return None
        return CompoundByShapes(parts).compound


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
        if assy is False:
            return None
        if isinstance(assy, Assembly):
            return assy
        try:
            return cls._all[assy]
        except KeyError:
            return cls.get_active()

    @classmethod
    def create_assy(cls, label, parent=None, active=True, *parts):
        """
        Create an assembly.

        :param label:
        :param parent:
        :param active:
        :param parts:

        :return:
        """
        if label in cls._all:
            return None
        if parent is None:
            parent = cls._master
        else:
            parent = cls.get_assy(parent)
        assy = Assembly(label, parent)
        if active:
            cls._active = assy
        assy.add_parts(*parts)
        cls._all[label] = assy
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
    def get_part(cls, label, assy=None):
        """
        Get a part from the assembly.
        
        :param label: 
        :param assy:
         
        :return: 
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return None
        return assy.get_part(label)

    @classmethod
    def get_parts(cls, assy=None, include_subassy=True):
        """
        Get parts from assembly.

        :param assy:
        :param include_subassy

        :return:
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return False
        return assy.get_parts(include_subassy)

    @classmethod
    def prepare_shape_to_mesh(cls, assy=None):
        """
        Prepare a shape for meshing.

        :param assy:

        :return:
        """
        if assy is None:
            assy = cls._master
        else:
            assy = cls.get_assy(assy)
        if not assy:
            return None
        return assy.prepare_shape_to_mesh()

    @classmethod
    def add_metadata(cls, key, value, assy=None):
        """
        Add metadata to the assembly.

        :param key:
        :param value:
        :param assy:

        :return:
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return False
        return assy.add_metadata(key, value)

    @classmethod
    def get_metadata(cls, key, assy=None):
        """
        Get metadata.

        :param key:
        :param assy:

        :return:
        """
        assy = cls.get_assy(assy)
        if not isinstance(assy, Assembly):
            return None
        return assy.get_metadata(key)
