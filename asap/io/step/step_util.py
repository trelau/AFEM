from OCC.IFSelect import IFSelect_RetError
from OCC.Interface import Interface_Static_SetCVal
from OCC.STEPControl import STEPControl_AsIs, STEPControl_Writer

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
