#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN AS "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR 
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

from OCCT.IFSelect import (IFSelect_ItemsByEntity, IFSelect_RetError,
                           IFSelect_RetDone)
from OCCT.Interface import Interface_Static
from OCCT.STEPControl import (STEPControl_AsIs, STEPControl_Reader,
                              STEPControl_Writer)

from afem.config import Settings, units_dict
from afem.topology.check import CheckShape

__all__ = ["StepExport", "StepImport"]


class StepExport(STEPControl_Writer):
    """
    Export shapes to a STEP file
    
    :param str schema: Define schema for STEP file ('AP203', or 'AP214').
    :param str units: Units to convert STEP file to.
    """

    def __init__(self, schema='AP203', units=None):
        super(StepExport, self).__init__()
        Interface_Static.SetCVal_('write.step.schema', schema)
        try:
            units = units_dict[units]
        except KeyError:
            units = Settings.units
        Interface_Static.SetCVal_('write.step.unit', units)

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
            status = self.Transfer(shape, STEPControl_AsIs)
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
        status = self.Write(fn)
        if int(status) < int(IFSelect_RetError):
            return True
        return False


class StepImport(STEPControl_Reader):
    """
    Import a STEP file.
    """

    def __init__(self):
        super(StepImport, self).__init__()
        self._shape = None

    @property
    def shape(self):
        """
        :return: The shape.
        """
        return self._shape

    def read(self, fn):
        """
        Read a STEP file.
        
        :param str fn: The full path to the file.
         
        :return: *True* if file was imported, *False* if not.
        :rtype: bool
        """
        # Read file.
        status = self.ReadFile(fn)
        if int(status) > int(IFSelect_RetDone):
            return False

        # Convert to desired units.
        Interface_Static.SetCVal_("xstep.cascade.unit", Settings.units)

        # Check
        self.PrintCheckLoad(False, IFSelect_ItemsByEntity)
        self.PrintCheckTransfer(False, IFSelect_ItemsByEntity)

        # Transfer
        self.TransferRoot(1)
        self._shape = self.Shape(1)
        return True
