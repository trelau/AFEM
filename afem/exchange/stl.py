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

from OCCT.StlAPI import StlAPI_Writer

__all__ = ["StlExport"]


class StlExport(StlAPI_Writer):
    """
    Export shape to STL file.
    """

    def __init__(self):
        super(StlExport, self).__init__()

    def write(self, shape, fn):
        """
        Converts shape to STL format and writes to a file.
        
        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.
        :param str fn: The filename.
         
        :return: None.
        """
        self.Write(shape, fn)
