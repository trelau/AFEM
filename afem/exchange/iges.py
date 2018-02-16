#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2018 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
from OCCT.IFSelect import IFSelect_RetDone
from OCCT.IGESControl import IGESControl_Reader, IGESControl_Writer
from OCCT.Interface import Interface_Static

from afem.config import Settings, units_dict

__all__ = ["IgesWrite", "IgesRead"]


class IgesWrite(object):
    """
    Write shape to an IGES file.

    :param str units: Units to convert IGES file to.
    :param int modecr: Option for writing faces. If 0, faces will be translated
        to IGES 144 (Trimmed Surface) entities. If 1, faces will be translated
        to IGES 510 (Face) entities and the IGES face will contain BRep
        entities.
    """

    def __init__(self, units='inch', modecr=0):
        self._writer = IGESControl_Writer(units, modecr)
        try:
            units = units_dict[units]
        except KeyError:
            units = Settings.units
        Interface_Static.SetCVal_('write.step.unit', units)

    @property
    def object(self):
        """
        :return: The IGES writer object.
        :rtype: OCCT.IGESControl.IGESControl_Writer
        """
        return self._writer

    def add_shape(self, shape):
        """
        Add the shape to the exported entities.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if shape was transferred, *False* if not.
        :rtype: bool
        """
        return self._writer.AddShape(shape)

    def add_geom(self, geom):
        """
        Add the geometry to the exported entities.

        :param OCCT.Geom.Geom_Geometry geom: The geometry.

        :return: *True* if shape was transferred, *False* if not.
        :rtype: bool
        """
        return self._writer.AddGeom(geom)

    def write(self, fn='afem.igs'):
        """
        Write the IGES file.

        :param str fn: The filename.

        :return: *True* if written, *False* if not.
        :rtype: bool
        """
        return self._writer.Write(fn)


class IgesRead(object):
    """
    Read an IGES file.

    :param str fn: The file to read.
    """

    def __init__(self, fn):
        self._reader = IGESControl_Reader()

        # Read file
        status = self._reader.ReadFile(fn)
        if status != IFSelect_RetDone:
            raise RuntimeError("Error reading IGES file.")

        # Convert to desired units
        Interface_Static.SetCVal_("xstep.cascade.unit", Settings.units)

        # Transfer
        nroots = self._reader.TransferRoots()
        if nroots > 0:
            self._shape = self._reader.OneShape()

    @property
    def object(self):
        """
        :return: The IGES reader object.
        :rtype: OCCT.IGESControl.IGESControl_Reader
        """
        return self._reader

    @property
    def shape(self):
        """
        :return: The main shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape
