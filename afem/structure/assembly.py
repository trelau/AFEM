from afem.structure.utils import order_parts_by_id
from afem.topology.create import CompoundByShapes

__all__ = ["Assembly", "AssemblyAPI"]


class Assembly(object):
    """
    Structural assembly. This is a minimal implementation only to help
    organize parts during the construction process.

    :param str label: The label.
    :param parent: The parent assembly, if any.
    :type parent: afem.structure.assembly.Assembly or None
    """

    def __init__(self, label, parent=None):
        self._label = label
        self._parent = parent
        self._children = set()
        self._parts = set()
        if isinstance(self._parent, Assembly):
            self._parent._children.add(self)
        self._metadata = {}

    @property
    def label(self):
        """
        :return: The label.
        :rtype: str
        """
        return self._label

    @property
    def parent(self):
        """
        :return: The parent assembly, if any.
        :rtype: afem.structure.assembly.Assembly or None
        """
        return self._parent

    @property
    def parts(self):
        """
        :return: List of all parts.
        :rtype: list[afem.structure.entities.Part]
        """
        return list(self._parts)

    @property
    def elements(self):
        """
        :return: The elements of all parts.
        :rtype: set(afem.fem.elements.Element)
        """
        elm_set = set()
        for part in self.parts:
            elm_set.update(part.elements)
        return elm_set

    @property
    def nodes(self):
        """
        :return: The nodes of all parts.
        :rtype: set(afem.fem.nodes.Node)
        """
        node_set = set()
        for part in self.parts:
            node_set.update(part.nodes)
        return node_set

    @property
    def metadata(self):
        """
        :return: The metadata dictionary.
        :rtype: dict
        """
        return self._metadata

    def add_metadata(self, key, value):
        """
        Add metadata to the assembly.

        :param key: The key.
        :param value: The value.

        :return: None.
        """
        self._metadata[key] = value

    def get_metadata(self, key):
        """
        Get metadata.

        :param key: The key.

        :return: The value.

        :raise KeyError: If the key is not in the dictionary.
        """
        return self._metadata[key]

    def activate(self):
        """
        Activate this assembly.

        :return: None
        """
        AssemblyAPI._active = self

    def add_parts(self, *parts):
        """
        Add parts to the assembly.

        :param afem.structure.entities.Part parts: Part(s) to add.

        :return: None.
        """
        part_set = set(parts)
        self._parts.update(part_set)

    def get_part(self, label):
        """
        Get a part in the assembly by label.

        :param str label: Part label.

        :return: The part.
        :rtype: afem.structure.entities.Part

        :raise KeyError: If the part is not found.
        """
        for part in self.parts:
            if part.label == label:
                return part
        raise KeyError('Part with given label could not be found in the '
                       'assembly.')

    def get_parts(self, include_subassy=True, rtype=None, order=False):
        """
        Get all the parts from the assembly and its sub-assemblies.

        :param bool include_subassy: Option to recursively include parts
            from any sub-assemblies.
        :param rtype: Option to return only parts of a certain type. Provide a
            class to check if the part is of the given type using
            *isinstance()*.
        :param bool order: Option to order parts by their ID.

        :return: List of parts.
        :rtype: list[afem.structure.entities.Part]
        """
        parts = []
        for part in self.parts:
            if rtype is None:
                parts.append(part)
            else:
                if isinstance(part, rtype):
                    parts.append(part)

        if include_subassy:
            for assy in self._children:
                parts += assy.get_parts(True, rtype)

        if not order:
            return parts
        return order_parts_by_id(parts)

    def prepare_shape_to_mesh(self, include_subassy=True):
        """
        Prepare a shape to mesh using the parts in the assembly and its
        sub-assemblies. This puts all the parts into a single compound which
        can be used as the master shape for the meshing process.

        :param bool include_subassy: Option to recursively include parts
            from any sub-assemblies.

        :return: The parts as a compound.
        :rtype: OCC.TopoDS.TopoDS_Compound
        """
        parts = self.get_parts(include_subassy)
        if not parts:
            return None
        return CompoundByShapes(parts).compound


class AssemblyAPI(object):
    """
    Assembly API. This stores all created assemblies so data can be accessed
    from one place. There is always a master model named '_master' that is
    created at initialization. This will be the parent of all assemblies
    that are not provided a parent when created. No assemblies should be
    labeled '_master'.
    """
    _master = Assembly('_master', None)
    _all = {'_master': _master}
    _active = _master

    @classmethod
    def get_master(cls):
        """
        Get the master assembly.

        :return: The master assembly.
        :rtype: afem.structure.assembly.Assembly
        """
        return cls._master

    @classmethod
    def get_active(cls):
        """
        Get the active assembly.

        :return: Active assembly.
        :rtype: afem.structure.assembly.Assembly
        """
        assy = cls._active
        if not isinstance(assy, Assembly):
            return cls._master
        return assy

    @classmethod
    def get_assy(cls, assy=None):
        """
        Get an assembly. If a string is provided then the assembly label
        will be used to get the assembly. If an assembly instance is
        provided then it is simply returned. If ``None`` is provided then
        the active assembly is returned.

        :param assy: Assembly to get.
        :type assy: str or afem.structure.assembly.Assembly or None

        :return: The assembly.
        :rtype: afem.structure.assembly.Assembly
        """
        if isinstance(assy, Assembly):
            return assy

        try:
            return cls._all[assy]
        except KeyError:
            return cls.get_active()

    @classmethod
    def make_active(cls, assy):
        """
        Activate the assembly. This method uses the *AssemblyAPI.get_assy()*
        method to get the assembly to activate.

        :param assy: Assembly to activate.
        :type assy: str or afem.structure.assembly.Assembly or None

        :return: None.
        """
        assy = cls.get_assy(assy)
        assy.activate()

    @classmethod
    def create_assy(cls, label, parent=None, active=True, *parts):
        """
        Create an assembly.

        :param str label: The label.
        :param parent: The parent assembly. If ``None`` then the active
            assembly is used.
        :type parent: str or afem.structure.assembly.Assembly or None
        :param bool active: Option to make the new assembly the active one.

        :param afem.structure.entities.Part parts: The part(s) to add to the
            assembly, if any.

        :return: New assembly.
        :rtype: afem.structure.assembly.Assembly
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

        :param assy: The assembly. If ``None`` then the active assembly is
            used.
        :type assy: str or afem.structure.assembly.Assembly or None
        :param afem.structure.entities.Part parts: The part(s) to add.

        :return: None.
        """
        assy = cls.get_assy(assy)
        assy.add_parts(*parts)

    @classmethod
    def get_part(cls, label, assy=None):
        """
        Get a part from the assembly using its label.
        
        :param str label: The part label.
        :param assy: The assembly. If ``None`` then the active assembly is
            used.
        :type assy: str or afem.structure.assembly.Assembly or None
         
        :return: The part.
        :rtype: afem.structure.entities.Part

        :raise KeyError: If the part is not found.
        """
        assy = cls.get_assy(assy)
        return assy.get_part(label)

    @classmethod
    def get_parts(cls, assy=None, include_subassy=True, rtype=None,
                  order=False):
        """
        Get parts from assembly.

        :param assy: The assembly. If ``None`` then the active assembly is
            used.
        :type assy: str or afem.structure.assembly.Assembly or None
        :param bool include_subassy: Option to recursively include parts
            from any sub-assemblies.
        :param rtype: Option to return only parts of a certain type. Provide a
            class to check if the part is of the given type using
            *isinstance()*.
        :param bool order: Option to order parts by their ID.

        :return: The parts.
        :rtype: list[afem.structure.entities.Part]
        """
        assy = cls.get_assy(assy)
        return assy.get_parts(include_subassy, rtype, order)

    @classmethod
    def prepare_shape_to_mesh(cls, assy='_master', include_subassy=True):
        """
        Prepare a shape to mesh using the parts in the assembly and its
        sub-assemblies. This puts all the parts into a single compound which
        can be used as the master shape for the meshing process.

        :param assy: The assembly. If ``None`` then the active assembly is
            used. By default the master model is used.
        :type assy: str or afem.structure.assembly.Assembly or None
        :param bool include_subassy: Option to recursively include parts
            from any sub-assemblies.

        :return: The parts as a compound.
        :rtype: OCC.TopoDS.TopoDS_Compound
        """
        assy = cls.get_assy(assy)
        return assy.prepare_shape_to_mesh(include_subassy)

    @classmethod
    def add_metadata(cls, key, value, assy=None):
        """
        Add metadata to the assembly.

        :param key: The key.
        :param value: The value.
        :param assy: The assembly. If ``None`` then the active assembly is
            used.
        :type assy: str or afem.structure.assembly.Assembly or None

        :return: None.
        """
        assy = cls.get_assy(assy)
        return assy.add_metadata(key, value)

    @classmethod
    def get_metadata(cls, key, assy=None):
        """
        Get metadata.

        :param key: They key.
        :param assy: The assembly. If ``None`` then the active assembly is
            used.
        :type assy: str or afem.structure.assembly.Assembly or None

        :return: The value.

        :raise KeyError: If the key is not in the dictionary.
        """
        assy = cls.get_assy(assy)
        return assy.get_metadata(key)
