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
from collections import defaultdict

from OCCT.STEPCAFControl import STEPCAFControl_Reader
from OCCT.TCollection import TCollection_ExtendedString
from OCCT.TDF import TDF_ChildIterator
from OCCT.TDataStd import TDataStd_Name
from OCCT.TDocStd import TDocStd_Document
from OCCT.TNaming import TNaming_NamedShape
from OCCT.XCAFApp import XCAFApp_Application
from OCCT.XCAFDoc import XCAFDoc_DocumentTool

from afem.topology.check import CheckShape
from afem.topology.create import CompoundByShapes


class XdeRead(object):
    """
    Read models using OpenCASCADE Extended Data Exchange.

    :param str fn: The full path to a file to read. If provided, the file
        extension will be used to guess which method to call. For example, if
        the file ended with 'stp' or 'step', then the *read_step()* method is
        called.

    :raise TypeError: If a file is provided but the extension is not recognized
        or supported.

    .. note::

        If multiple shapes are found with the same label then they are placed
        into a single TopoDS_Compound.
    """

    def __init__(self, fn=None):
        fmt = TCollection_ExtendedString('MDTV-Standard')
        self._doc = TDocStd_Document(fmt)
        self._app = XCAFApp_Application.GetApplication_()
        self._app.InitDocument(self._doc)

        self._labels_to_shapes = {}
        self._shape = None

        if fn is not None:
            if fn.endswith('.step') or fn.endswith('.stp'):
                self.read_step(fn)
            else:
                raise TypeError('File extension not supported.')

    @property
    def shape(self):
        """
        :return: The main shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape

    @property
    def labels(self):
        """
        :return: List of labels.
        :rtype: list(str)
        """
        return list(self._labels_to_shapes.keys())

    @property
    def label_shape_iter(self):
        """
        :return: Yield pairs of (label, shape).
        :rtype: tuple(str, OCCT.TopoDS.TopoDS_Shape)
        """
        for label in self._labels_to_shapes:
            shape = self._labels_to_shapes[label]
            yield label, shape

    def _transfer(self):
        """
        Transfer named shapes.
        """
        label = XCAFDoc_DocumentTool.ShapesLabel_(self._doc.Main())
        iter_ = TDF_ChildIterator(label, True)
        names_to_shapes = defaultdict(list)
        while iter_.More():
            current = iter_.Value()
            iter_.Next()

            name, shape = None, None
            # Name
            attr1 = TDataStd_Name()
            status1, attr1 = current.FindAttribute(attr1.GetID_(), attr1)
            if status1:
                name = attr1.Get().ToExtString()

            # Named shape
            attr2 = TNaming_NamedShape()
            status2, attr2 = current.FindAttribute(attr2.GetID_(), attr2)
            if status2:
                shape = CheckShape.to_shape(attr2.Get())

            if not status1 or not status2:
                continue

            names_to_shapes[name].append(shape)

        for name in names_to_shapes:
            shapes = names_to_shapes[name]
            if len(shapes) > 1:
                shape = CompoundByShapes(shapes).compound
            else:
                shape = shapes[0]
            self._labels_to_shapes[name] = shape

        return True

    def read_step(self, fn):
        """
        Read a STEP file.

        :param str fn: The filename.

        :return: *True* if successful, *False* if not.

        :raise RuntimeError: If the file cannot be read.
        """
        reader = STEPCAFControl_Reader()
        reader.SetNameMode(True)
        status = reader.Perform(fn, self._doc)
        if not status:
            raise RuntimeError("Error reading STEP file.")

        self._shape = reader.Reader().OneShape()
        return self._transfer()

    def get_shape(self, label):
        """
        Get a shape by its label.

        :param str label: Shape label.

        :return: The shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape

        :raise KeyError: If ``label`` is not present.
        """
        return self._labels_to_shapes[label]
