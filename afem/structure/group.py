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
from OCCT.BinXCAFDrivers import BinXCAFDrivers
from OCCT.PCDM import PCDM_StoreStatus
from OCCT.TCollection import TCollection_AsciiString, TCollection_ExtendedString
from OCCT.TDF import TDF_ChildIterator
from OCCT.TDataStd import TDataStd_Name, TDataStd_AsciiString
from OCCT.TDocStd import TDocStd_Document
from OCCT.TNaming import TNaming_NamedShape
from OCCT.XCAFApp import XCAFApp_Application
from OCCT.XCAFDoc import XCAFDoc_DocumentTool, XCAFDoc_Color

from afem.structure.utils import order_parts_by_id
from afem.topology.check import CheckShape
from afem.topology.create import CompoundByShapes, EdgeByCurve, FaceBySurface
from afem.topology.explore import ExploreShape

__all__ = ["Group", "GroupAPI"]


class Group(object):
    """
    Group of parts.

    :param str label: The label.
    :param parent: The parent group, if any.
    :type parent: afem.structure.group.Group or None
    """

    def __init__(self, label, parent=None):
        self._label = label
        self._parent = parent
        self._children = set()
        self._parts = set()
        if isinstance(self._parent, Group):
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
        :return: The parent group, if any.
        :rtype: afem.structure.group.Group or None
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
    def metadata(self):
        """
        :return: The metadata dictionary.
        :rtype: dict
        """
        return self._metadata

    def add_metadata(self, key, value):
        """
        Add metadata to the group.

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
        Activate this group.

        :return: None
        """
        GroupAPI._active = self

    def add_parts(self, *parts):
        """
        Add parts to the group.

        :param afem.structure.entities.Part parts: Part(s) to add.

        :return: None.
        """
        part_set = set(parts)
        self._parts.update(part_set)

    def get_part(self, label):
        """
        Get a part in the group by label.

        :param str label: Part label.

        :return: The part.
        :rtype: afem.structure.entities.Part

        :raise KeyError: If the part is not found.
        """
        for part in self.parts:
            if part.label == label:
                return part
        raise KeyError('Part with given label could not be found in the '
                       'group.')

    def get_parts(self, include_subgroup=True, rtype=None, order=False):
        """
        Get all the parts from the group and its subgroups.

        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.
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

        if include_subgroup:
            for group in self._children:
                parts += group.get_parts(True, rtype)

        if not order:
            return parts
        return order_parts_by_id(parts)

    def prepare_shape_to_mesh(self, include_subgroup=True):
        """
        Prepare a shape to mesh using the parts in the group and its
        subgroups. This puts all the parts into a single compound which
        can be used as the master shape for the meshing process.

        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.

        :return: The parts as a compound.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        return self.as_compound(include_subgroup)

    def as_compound(self, include_subgroup=True):
        """
        Build a TopoDS_Compound from all the parts of the group.

        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.

        :return: The parts as a compound.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        parts = self.get_parts(include_subgroup)
        shapes = [part.shape for part in parts]
        return CompoundByShapes(shapes).compound


class GroupAPI(object):
    """
    Group API. This stores all created groups so data can be accessed
    from one place. There is always a master model named '_master' that is
    created at initialization. This will be the parent of all groups
    that are not provided a parent when created. No groups should be
    labeled '_master'.
    """
    _master = Group('_master', None)
    _all = {'_master': _master}
    _active = _master

    @classmethod
    def reset(cls):
        """
        Reset master group and data structure and reset Part index back to 1.
        This should delete all groups unless they are referenced somewhere
        else.

        :return: None.
        """
        cls._master = Group('_master', None)
        cls._all = {'_master': cls._master}
        cls._active = cls._master

        from afem.structure.entities import Part

        Part.reset()

    @classmethod
    def get_master(cls):
        """
        Get the master group.

        :return: The master group.
        :rtype: afem.structure.group.Group
        """
        return cls._master

    @classmethod
    def get_active(cls):
        """
        Get the active group.

        :return: Active group.
        :rtype: afem.structure.group.Group
        """
        group = cls._active
        if not isinstance(group, Group):
            return cls._master
        return group

    @classmethod
    def get_group(cls, group=None):
        """
        Get an group. If a string is provided then the group label
        will be used to get the group. If a group instance is
        provided then it is simply returned. If ``None`` is provided then
        the active group is returned.

        :param group: Group to get.
        :type group: str or afem.structure.group.Group or None

        :return: The group.
        :rtype: afem.structure.group.Group
        """
        if isinstance(group, Group):
            return group

        try:
            return cls._all[group]
        except KeyError:
            return cls.get_active()

    @classmethod
    def make_active(cls, group):
        """
        Activate the group. This method uses the *GroupAPI.get_group()*
        method to get the group to activate.

        :param group: Group to activate.
        :type group: str or afem.structure.group.Group or None

        :return: None.
        """
        group = cls.get_group(group)
        group.activate()

    @classmethod
    def create_group(cls, label, parent=None, active=True, *parts):
        """
        Create an group.

        :param str label: The label.
        :param parent: The parent group. If ``None`` then the master
            group is used.
        :type parent: str or afem.structure.group.Group or None
        :param bool active: Option to make the new group the active one.

        :param afem.structure.entities.Part parts: The part(s) to add to the
            group, if any.

        :return: New group.
        :rtype: afem.structure.group.Group
        """
        if label in cls._all:
            return None
        if parent is None:
            parent = cls._master
        else:
            parent = cls.get_group(parent)
        group = Group(label, parent)
        if active:
            cls._active = group
        group.add_parts(*parts)
        cls._all[label] = group
        return group

    @classmethod
    def add_parts(cls, group, *parts):
        """
        Add parts to the group.

        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None
        :param afem.structure.entities.Part parts: The part(s) to add.

        :return: None.
        """
        group = cls.get_group(group)
        group.add_parts(*parts)

    @classmethod
    def get_part(cls, label, group=None):
        """
        Get a part from the group using its label.
        
        :param str label: The part label.
        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None
         
        :return: The part.
        :rtype: afem.structure.entities.Part

        :raise KeyError: If the part is not found.
        """
        group = cls.get_group(group)
        return group.get_part(label)

    @classmethod
    def get_parts(cls, group=None, include_subgroup=True, rtype=None,
                  order=False):
        """
        Get parts from group.

        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None
        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.
        :param rtype: Option to return only parts of a certain type. Provide a
            class to check if the part is of the given type using
            *isinstance()*.
        :param bool order: Option to order parts by their ID.

        :return: The parts.
        :rtype: list[afem.structure.entities.Part]
        """
        group = cls.get_group(group)
        return group.get_parts(include_subgroup, rtype, order)

    @classmethod
    def prepare_shape_to_mesh(cls, group='_master', include_subgroup=True):
        """
        Prepare a shape to mesh using the parts in the group and its
        subgroups. This puts all the parts into a single compound which
        can be used as the master shape for the meshing process.

        :param group: The group. If ``None`` then the active group is
            used. By default the master model is used.
        :type group: str or afem.structure.group.Group or None
        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.

        :return: The parts as a compound.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        group = cls.get_group(group)
        return group.prepare_shape_to_mesh(include_subgroup)

    @classmethod
    def as_compound(cls, group='_master', include_subgroup=True):
        """
        Build a TopoDS_Compound from all the parts of the group.

        :param group: The group. If ``None`` then the active group is
            used. By default the master model is used.
        :type group: str or afem.structure.group.Group or None
        :param bool include_subgroup: Option to recursively include parts
            from any subgroups.

        :return: The parts as a compound.
        :rtype: OCCT.TopoDS.TopoDS_Compound
        """
        group = cls.get_group(group)
        return group.as_compound(include_subgroup)

    @classmethod
    def add_metadata(cls, key, value, group=None):
        """
        Add metadata to the group.

        :param key: The key.
        :param value: The value.
        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None

        :return: None.
        """
        group = cls.get_group(group)
        return group.add_metadata(key, value)

    @classmethod
    def get_metadata(cls, key, group=None):
        """
        Get metadata.

        :param key: They key.
        :param group: The group. If ``None`` then the active group is
            used.
        :type group: str or afem.structure.group.Group or None

        :return: The value.

        :raise KeyError: If the key is not in the dictionary.
        """
        group = cls.get_group(group)
        return group.get_metadata(key)

    @classmethod
    def save_model(cls, fn):
        """
        Save the model.

        :param str fn: The filename. The extension will be ".xbf" and appended
            if not provided.

        :return: *True* if saved, *False* otherwise.
        :rtype: bool
        """
        group = cls.get_master()

        # Create document and application
        fmt = TCollection_ExtendedString('BinXCAF')
        doc = TDocStd_Document(fmt)
        app = XCAFApp_Application.GetApplication_()
        BinXCAFDrivers.DefineFormat_(app)
        app.InitDocument(doc)

        # Store parts as top-level shapes
        # TODO Support group hierarchy
        tool = XCAFDoc_DocumentTool.ShapeTool_(doc.Main())
        parts = group.get_parts()
        for part in parts:
            # New label
            label = tool.AddShape(part.shape, False)

            # Name
            txt = TCollection_ExtendedString(part.label)
            TDataStd_Name.Set_(label, txt)

            # Type
            type_ = TCollection_AsciiString(part.type)
            TDataStd_AsciiString.Set_(label, type_)

            # Color
            XCAFDoc_Color.Set_(label, part.color)

            # Reference curve
            if part.has_cref:
                edge = EdgeByCurve(part.cref).edge
                label = tool.AddShape(edge, False)
                txt = TCollection_ExtendedString(part.label)
                TDataStd_Name.Set_(label, txt)
                type_ = TCollection_AsciiString('CREF')
                TDataStd_AsciiString.Set_(label, type_)

            # Reference surface
            if part.has_sref:
                face = FaceBySurface(part.sref).face
                label = tool.AddShape(face, False)
                txt = TCollection_ExtendedString(part.label)
                TDataStd_Name.Set_(label, txt)
                type_ = TCollection_AsciiString('SREF')
                TDataStd_AsciiString.Set_(label, type_)

        if not fn.endswith('.xbf'):
            fn += '.xbf'
        txt = TCollection_ExtendedString(fn)
        status = app.SaveAs(doc, txt)
        return status == PCDM_StoreStatus.PCDM_SS_OK

    @classmethod
    def load_model(cls, fn, group=None):
        """
        Load a model.

        :param str fn: The filename. The extension will be ".xbf" and appended
            if not provided.
        :param afem.structure.group.Group group: The group to load the parts
            into. If *None* then the active group is used.

        :return: *True* if loaded, *False* otherwise.
        :rtype: bool
        """
        from afem.structure.create import CreatePartByName

        group = cls.get_group(group)

        # Create document and application
        fmt = TCollection_ExtendedString('BinXCAF')
        doc = TDocStd_Document(fmt)
        app = XCAFApp_Application.GetApplication_()
        BinXCAFDrivers.DefineFormat_(app)
        app.InitDocument(doc)

        # Open the document
        if not fn.endswith('.xbf'):
            fn += '.xbf'
        txt = TCollection_ExtendedString(fn)
        status, doc = app.Open(txt, doc)

        # Get the main label and iterate on top-level children which
        # should be parts
        label = XCAFDoc_DocumentTool.ShapesLabel_(doc.Main())
        iter_ = TDF_ChildIterator(label)
        part_data = []
        cref_to_part = {}
        sref_to_part = {}
        while iter_.More():
            current = iter_.Value()
            iter_.Next()

            name, shape, type_, color = 4 * [None]

            # Name
            attr1 = TDataStd_Name()
            status1, attr1 = current.FindAttribute(attr1.GetID_(), attr1)
            if status1:
                name = attr1.Get().ToExtString()

            # Named shape
            attr2 = TNaming_NamedShape()
            status2, attr2 = current.FindAttribute(attr2.GetID_(), attr2)
            if status2:
                shape = attr2.Get()

            # Type
            attr3 = TDataStd_AsciiString()
            status3, attr3 = current.FindAttribute(attr3.GetID_(), attr3)
            if status3:
                type_ = attr3.Get().ToCString()

            # Color
            attr4 = XCAFDoc_Color()
            status4, attr4 = current.FindAttribute(attr4.GetID_(), attr4)
            if status4:
                color = attr4.GetColor()

            if None in [name, shape, type_]:
                continue

            # Check for reference geometry
            if type_ == 'CREF':
                edge = CheckShape.to_edge(shape)
                cref_to_part[name] = ExploreShape.curve_of_edge(edge)
                continue

            if type_ == 'SREF':
                face = CheckShape.to_face(shape)
                sref_to_part[name] = ExploreShape.surface_of_face(face)
                continue

            # Add part data
            part_data.append((type_, name, shape, color))

        # Create parts
        # TODO Support group hierarchy
        for type_, label, shape, color in part_data:
            cref, sref = None, None
            if label in cref_to_part:
                cref = cref_to_part[label]
            if label in sref_to_part:
                sref = sref_to_part[label]

            part = CreatePartByName(type_, label=label, shape=shape,
                                    group=group, cref=cref, sref=sref).part
            if color is not None:
                part.color = color

        return True
