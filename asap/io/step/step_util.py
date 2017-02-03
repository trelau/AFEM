from OCC.IFSelect import IFSelect_RetError
from OCC.Interface import Interface_Static_SetCVal
from OCC.STEPControl import STEPControl_AsIs, STEPControl_Writer

from ...topology import ShapeTools


class StepExport(STEPControl_Writer):
    """
    Export shapes to a step file
    """

    def __init__(self, schema='AP203'):
        super(StepExport, self).__init__()
        Interface_Static_SetCVal('write.step.schema', schema)

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
