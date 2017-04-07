from OCC.IFSelect import IFSelect_ItemsByEntity, IFSelect_RetError
from OCC.Interface import Interface_Static_SetCVal
from OCC.STEPControl import STEPControl_AsIs, STEPControl_Reader, \
    STEPControl_Writer

from ...config import Settings, units_dict
from ...topology import ShapeTools


class StepExport(STEPControl_Writer):
    """
    Export shapes to a STEP file
    
    :param str schema: Define schema for STEP file ('AP203', or 'AP214').
    :param str units: Units to convert STEP file to.
    """

    def __init__(self, schema='AP203', units=None):
        super(StepExport, self).__init__()
        Interface_Static_SetCVal('write.step.schema', schema)
        try:
            units = units_dict[units]
        except KeyError:
            units = Settings.units
        Interface_Static_SetCVal('write.step.unit', units)

    def transfer(self, *shapes):
        """
        Transfer and add the shapes to the exported entities.

        :param shapes:

        :return:
        """
        added_shape = False
        for shape in shapes:
            shape = ShapeTools.to_shape(shape)
            if not shape:
                continue
            status = self.Transfer(shape, STEPControl_AsIs)
            if status < IFSelect_RetError:
                added_shape = True
        return added_shape

    def write(self, fn='asap.stp'):
        """
        Write STEP file.

        :param fn:

        :return:
        """
        status = self.Write(fn)
        if status < IFSelect_RetError:
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
        return self._shape

    def read(self, fn):
        """
        Read a STEP file.
        
        :param fn:
         
        :return: 
        """
        # Read file.
        status = self.ReadFile(fn)
        if status > 1:
            return False

        # Convert to desired units.
        Interface_Static_SetCVal("xstep.cascade.unit", Settings.units)

        # Check
        self.PrintCheckLoad(False, IFSelect_ItemsByEntity)
        self.PrintCheckTransfer(False, IFSelect_ItemsByEntity)

        # Transfer
        self.TransferRoot(1)
        self._shape = self.Shape(1)
        return True
