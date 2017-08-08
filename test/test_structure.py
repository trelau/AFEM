import unittest

from afem.geometry import *
from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.structure.create import *
from afem.structure.entities import *
from afem.topology import *


class TestStructure(unittest.TestCase):
    ImportVSP.step_file('./test_io/vsp_transport_wing_split_350.stp')

    wing1 = ImportVSP.get_body('Wing')
    wing1.set_transparency(0.5)

    def _show(self, *items):
        for item in items:
            Viewer.add(item)
        Viewer.show()

    def test_spar_by_parameters(self):
        builder = SparByParameters('spar', 0.5, 0., 0.5, 0.75, self.wing1)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_spar_by_points(self):
        p1 = self.wing1.eval(0.5, 0.)
        p2 = self.wing1.eval(0.5, 0.75)
        builder = SparByPoints('spar', p1, p2, self.wing1)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_spar_by_surface(self):
        p0 = self.wing1.eval(0.5, 0.)
        pln = PlaneByAxes(p0, 'yz').plane
        builder = SparBySurface('spar', pln, self.wing1)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_spar_between_shapes(self):
        p = self.wing1.eval(0.5, 0.)
        pln1 = PlaneByAxes(p, 'xz').plane
        shape1 = FaceBySurface(pln1).face
        p = self.wing1.eval(0.5, 0.5)
        pln2 = PlaneByAxes(p, 'xz').plane
        shape2 = FaceBySurface(pln2).face
        p = self.wing1.eval(0.5, 0.)
        pln3 = PlaneByNormal(p, (0.5, -0.2, 0.)).plane
        basis_shape = FaceBySurface(pln3).face
        builder = SparBetweenShapes('spar', shape1, shape2, self.wing1,
                                    basis_shape)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_rib_by_parameters(self):
        builder = RibByParameters('rib', 0.15, 0.15, 0.65, 0.15, self.wing1)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_by_points(self):
        p1 = self.wing1.eval(0.15, 0.15)
        p2 = self.wing1.eval(0.65, 0.15)
        builder = RibByPoints('rib', p1, p2, self.wing1)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_by_surface(self):
        p0 = self.wing1.eval(0.5, 0.5)
        pln = PlaneByAxes(p0, 'xz').plane
        builder = RibBySurface('rib', pln, self.wing1)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_between_shapes(self):
        p = self.wing1.eval(0.15, 0.5)
        pln1 = PlaneByAxes(p, 'yz').plane
        shape1 = FaceBySurface(pln1).face
        p = self.wing1.eval(0.65, 0.5)
        pln2 = PlaneByAxes(p, 'yz').plane
        shape2 = FaceBySurface(pln2).face
        p = self.wing1.eval(0.5, 0.5)
        pln3 = PlaneByAxes(p, 'xz').plane
        basis_shape = FaceBySurface(pln3).face
        builder = RibBetweenShapes('rib', shape1, shape2, self.wing1,
                                   basis_shape)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_between_planes_by_number(self):
        builder = SparByParameters('fspar', 0.15, 0.15, 0.15, 0.5, self.wing1)
        fspar = builder.spar
        builder = SparByParameters('rspar', 0.65, 0.15, 0.65, 0.5, self.wing1)
        rspar = builder.spar
        pln1 = PlaneByAxes(fspar.p1, 'xz').plane
        pln2 = PlaneByAxes(fspar.p2, 'xz').plane
        builder = RibsBetweenPlanesByNumber('rib', pln1, pln2, 10, fspar,
                                            rspar, self.wing1)
        self.assertEqual(builder.nribs, 10)
        self.assertEqual(builder.next_index, 11)
        for rib in builder.ribs:
            self.assertIsInstance(rib, Rib)

    def test_rib_between_planes_by_distance(self):
        builder = SparByParameters('fspar', 0.15, 0.15, 0.15, 0.5, self.wing1)
        fspar = builder.spar
        builder = SparByParameters('rspar', 0.65, 0.15, 0.65, 0.5, self.wing1)
        rspar = builder.spar
        pln1 = PlaneByAxes(fspar.p1, 'xz').plane
        pln2 = PlaneByAxes(fspar.p2, 'xz').plane
        builder = RibsBetweenPlanesByDistance('rib', pln1, pln2, 36., fspar,
                                              rspar, self.wing1)
        self.assertEqual(builder.nribs, 12)
        self.assertAlmostEqual(builder.spacing, 35.0892460667128, delta=0.001)
        self.assertEqual(builder.next_index, 13)
        for rib in builder.ribs:
            self.assertIsInstance(rib, Rib)


if __name__ == '__main__':
    unittest.main()
