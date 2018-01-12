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

from OCCT.IFSelect import (IFSelect_RetError,
                           IFSelect_RetDone)
from OCCT.Interface import Interface_Static
from OCCT.STEPControl import (STEPControl_AsIs, STEPControl_Writer,
                              STEPControl_Reader)

from afem.config import Settings, units_dict
from afem.topology.check import CheckShape

__all__ = ["StepWrite", "StepRead"]


class StepWrite(object):
    """
    Write shape to a STEP file.
    
    :param str schema: Define schema for STEP file ('AP203', or 'AP214').
    :param str units: Units to convert STEP file to.
    """

    def __init__(self, schema='AP203', units=None):
        self._writer = STEPControl_Writer()
        Interface_Static.SetCVal_('write.step.schema', schema)
        try:
            units = units_dict[units]
        except KeyError:
            units = Settings.units
        Interface_Static.SetCVal_('write.step.unit', units)

    @property
    def object(self):
        """
        :return: The STEP writer object.
        :rtype: OCCT.STEPControl.STEPControl_Writer
        """
        return self._writer

    def transfer(self, *shapes):
        """
        Transfer and add the shapes to the exported entities.

        :param OCCT.TopoDS.TopoDS_Shape shapes: The shape(s).

        :return: *True* if shape was transferred, *False* if not.
        :rtype: bool
        """
        added_shape = False
        for shape in shapes:
            shape = CheckShape.to_shape(shape)
            if not shape:
                continue
            status = self._writer.Transfer(shape, STEPControl_AsIs)
            if int(status) < int(IFSelect_RetError):
                added_shape = True
        return added_shape

    def write(self, fn='afem.stp'):
        """
        Write the STEP file.

        :param str fn: The filename.

        :return: *True* if written, *False* if not.
        :rtype: bool
        """
        status = self._writer.Write(fn)
        if int(status) < int(IFSelect_RetError):
            return True
        return False


class StepRead(object):
    """
    Read a STEP file.

    :param str fn: The file to read.
    """

    def __init__(self, fn):
        self._reader = STEPControl_Reader()

        # Read file
        status = self._reader.ReadFile(fn)
        if status != IFSelect_RetDone:
            raise RuntimeError("Error reading STEP file.")

        # Transfer
        nroots = self._reader.TransferRoots()
        if nroots > 0:
            self._shape = self._reader.OneShape()

    @property
    def object(self):
        """
        :return: The STEP reader object.
        :rtype: OCCT.STEPControl.STEPControl_Reader
        """
        return self._reader

    @property
    def shape(self):
        """
        :return: The main shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._shape
