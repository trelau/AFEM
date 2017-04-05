from OCC.BRepBndLib import brepbndlib_Add
from OCC.Bnd import Bnd_Box

from .checker import CheckGeom


class BBox(Bnd_Box):
    """
    Bounding box in 3-D.
    """

    def __init__(self):
        super(BBox, self).__init__()

    @property
    def is_void(self):
        return self.IsVoid()

    @property
    def pmin(self):
        if self.is_void:
            return None
        return CheckGeom.to_point(self.CornerMin())

    @property
    def pmax(self):
        if self.is_void:
            return None
        return CheckGeom.to_point(self.CornerMax())

    @property
    def xmin(self):
        if self.is_void:
            return None
        return self.CornerMin().X()

    @property
    def xmax(self):
        if self.is_void:
            return None
        return self.CornerMax().X()

    @property
    def ymin(self):
        if self.is_void:
            return None
        return self.CornerMin().Y()

    @property
    def ymax(self):
        if self.is_void:
            return None
        return self.CornerMax().Y()

    @property
    def zmin(self):
        if self.is_void:
            return None
        return self.CornerMin().Z()

    @property
    def zmax(self):
        if self.is_void:
            return None
        return self.CornerMax().Z()

    def add_point(self, pnt):
        """
        Add the point to the box.
        
        :param pnt:
         
        :return: 
        """
        pnt = CheckGeom.to_point(pnt)
        if not pnt:
            return False
        self.Add(pnt)
        return True

    def add_box(self, box):
        """
        Add the other box to the box.

        :param box:

        :return: 
        """
        if not isinstance(box, Bnd_Box):
            return False
        self.Add(box)
        return True

    def add_shape(self, shape):
        """
        Add shape to the boudning box.
        
        :param shape:
        
        :return: 
        """
        brepbndlib_Add(shape, self, True)
        return True

    def distance(self, box):
        """
        Calculate distance to other box.
        
        :param box:
         
        :return: 
        """
        if self.is_void or box.IsVoid():
            return None
        return self.Distance(box)
