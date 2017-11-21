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
