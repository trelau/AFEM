from .part import Part
from ..topology import ShapeTools


class CurvePart(Part):
    """
    Base class for curve-based parts.
    """

    def __init__(self, name, curve_shape):
        super(CurvePart, self).__init__(name)
        curve_shape = ShapeTools.to_shape(curve_shape)
        self.set_shape(curve_shape)
        self._cref = None
        self._set_cref()

    @property
    def edges(self):
        return ShapeTools.get_edges(self)

    @property
    def nedges(self):
        return len(self.edges)

    @property
    def cref(self):
        return self._cref

    @property
    def reshapes(self):
        return self.edges

    def _set_cref(self):
        """
        Set part reference curve if available. 
        """
        if self.IsNull():
            return False

        # Use longest wire to set cref.
        wires = ShapeTools.get_wires(self)
        if not wires:
            edges = ShapeTools.get_edges(self)
            wires = ShapeTools.connect_edges(edges)
        if not wires:
            return False
        wire = ShapeTools.longest_wire(wires)
        if not wire:
            return False

        # Concatenate to a single (possibly C0) curve.
        edge = ShapeTools.concatenate_wire(wire)
        if not edge:
            return False
        cref = ShapeTools.curve_of_edge(edge)
        if cref:
            self._cref = cref
            return True
        return False
