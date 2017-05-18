from OCC.TopoDS import TopoDS_Solid

from ..graphics.viewer import ViewableItem
from ..topology import ShapeTools


class Body(TopoDS_Solid, ViewableItem):
    """
    Base class for OML bodies.
    """

    def __init__(self, shape, name=None):
        super(Body, self).__init__()
        ViewableItem.__init__(self)
        self.set_shape(shape)
        self._metadata = {}
        self._name = name

    @property
    def shape(self):
        return self

    @shape.setter
    def shape(self, shape):
        self.set_shape(shape)

    @property
    def shell(self):
        return ShapeTools.outer_shell(self)

    @property
    def metadata(self):
        return self._metadata

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self.set_name(name)

    def set_name(self, name):
        """
        Set name of body.

        :param name:

        :return:
        """
        self._name = name

    def add_metadata(self, key, value):
        """
        Add metadata to the body.

        :param key:
        :param value:

        :return:
        """
        self._metadata[key] = value

    def get_metadata(self, key):
        """
        Get metadata.

        :param key:

        :return:
        """
        try:
            return self._metadata[key]
        except KeyError:
            return None

    def set_shape(self, shape):
        """
        Set the shape for the OML body.

        :param shape:

        :return:
        """
        shape = ShapeTools.to_solid(shape)
        if not shape:
            return False

        self.TShape(shape.TShape())
        self.Location(shape.Location())
        self.Orientation(shape.Orientation())
        return True

    def _bop(self, body, inplace, bop):
        """
        Perform the BOP.
        """
        if inplace:
            failed = False
        else:
            failed = None
        if bop == 'fuse':
            solids = ShapeTools.bfuse(self, body, 'solid')
        elif bop == 'common':
            solids = ShapeTools.bcommon(self, body, 'solid')
        elif bop == 'cut':
            solids = ShapeTools.bcut(self, body, 'solid')
        else:
            return failed
        if not solids:
            return failed
        solid = solids[0]
        if inplace:
            return self.set_shape(solid)
        return Body(solid)

    def fuse(self, body, inplace=False):
        """
        Fuse the bodies.

        :param body:
        :param inplace:

        :return:
        """
        return self._bop(body, inplace, 'fuse')

    def common(self, body, inplace=False):
        """
        Find the common solid between the bodies.

        :param body:
        :param inplace:

        :return:
        """
        return self._bop(body, inplace, 'common')

    def cut(self, body, inplace=False):
        """
        Cut this body with the other.

        :param body:
        :param inplace:

        :return:
        """
        return self._bop(body, inplace, 'cut')

    def section(self, body, rtype=''):
        """
        Find intersection between the bodies.

        :param body:
        :param rtype:

        :return:
        """
        return ShapeTools.bsection(self, body, rtype)

    def bop_algo(self, bodies, operation='fuse'):
        """
        Perform BOP on multiple bodies.

        :param operation:
        :param bodies:

        :return:
        """
        return ShapeTools.bop_algo([self], bodies, operation)
