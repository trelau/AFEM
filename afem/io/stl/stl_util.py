from OCC.StlAPI import StlAPI_Writer

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
        
        :param shape: 
        :param fn:
         
        :return: 
        """
        self.Write(shape, fn)
