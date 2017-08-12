import unittest

from afem.geometry import *
from afem.graphics import Viewer
from afem.io import ImportVSP
from afem.structure.create import *
from afem.structure.entities import *
from afem.topology import *


def _show(*items):
    for item in items:
        Viewer.add(item)
    Viewer.show()


class TestStructure(unittest.TestCase):
    ImportVSP.step_file('./test_io/777-200LR.stp')

    wing = ImportVSP.get_body('Wing')
    fuselage = ImportVSP.get_body('Fuselage')
    wing.set_transparency(0.5)
    fuselage.set_transparency(0.5)

    def test_spar_by_parameters(self):
        builder = SparByParameters('spar', 0.5, 0., 0.5, 0.75, self.wing)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_spar_by_points(self):
        p1 = self.wing.eval(0.5, 0.)
        p2 = self.wing.eval(0.5, 0.75)
        builder = SparByPoints('spar', p1, p2, self.wing)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_spar_by_surface(self):
        p0 = self.wing.eval(0.5, 0.)
        pln = PlaneByAxes(p0, 'yz').plane
        builder = SparBySurface('spar', pln, self.wing)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_spar_between_shapes(self):
        p = self.wing.eval(0.5, 0.)
        pln1 = PlaneByAxes(p, 'xz').plane
        shape1 = FaceBySurface(pln1).face
        p = self.wing.eval(0.5, 0.5)
        pln2 = PlaneByAxes(p, 'xz').plane
        shape2 = FaceBySurface(pln2).face
        p = self.wing.eval(0.5, 0.)
        pln3 = PlaneByNormal(p, (0.5, -0.2, 0.)).plane
        basis_shape = FaceBySurface(pln3).face
        builder = SparBetweenShapes('spar', shape1, shape2, self.wing,
                                    basis_shape)
        spar = builder.spar
        self.assertIsInstance(spar, Spar)

    def test_rib_by_parameters(self):
        builder = RibByParameters('rib', 0.15, 0.15, 0.65, 0.15, self.wing)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_by_points(self):
        p1 = self.wing.eval(0.15, 0.15)
        p2 = self.wing.eval(0.65, 0.15)
        builder = RibByPoints('rib', p1, p2, self.wing)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_by_surface(self):
        p0 = self.wing.eval(0.5, 0.5)
        pln = PlaneByAxes(p0, 'xz').plane
        builder = RibBySurface('rib', pln, self.wing)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_rib_between_shapes(self):
        p = self.wing.eval(0.15, 0.5)
        pln1 = PlaneByAxes(p, 'yz').plane
        shape1 = FaceBySurface(pln1).face
        p = self.wing.eval(0.65, 0.5)
        pln2 = PlaneByAxes(p, 'yz').plane
        shape2 = FaceBySurface(pln2).face
        p = self.wing.eval(0.5, 0.5)
        pln3 = PlaneByAxes(p, 'xz').plane
        basis_shape = FaceBySurface(pln3).face
        builder = RibBetweenShapes('rib', shape1, shape2, self.wing,
                                   basis_shape)
        rib = builder.rib
        self.assertIsInstance(rib, Rib)

    def test_ribs_between_planes_by_number(self):
        builder = SparByParameters('fspar', 0.15, 0.15, 0.15, 0.5, self.wing)
        fspar = builder.spar
        builder = SparByParameters('rspar', 0.65, 0.15, 0.65, 0.5, self.wing)
        rspar = builder.spar
        pln1 = PlaneByAxes(fspar.p1, 'xz').plane
        pln2 = PlaneByAxes(fspar.p2, 'xz').plane
        builder = RibsBetweenPlanesByNumber('rib', pln1, pln2, 10, fspar,
                                            rspar, self.wing)
        self.assertEqual(builder.nribs, 10)
        self.assertEqual(builder.next_index, 11)
        for rib in builder.ribs:
            self.assertIsInstance(rib, Rib)

    def test_ribs_between_planes_by_distance(self):
        builder = SparByParameters('fspar', 0.15, 0.15, 0.15, 0.5, self.wing)
        fspar = builder.spar
        builder = SparByParameters('rspar', 0.65, 0.15, 0.65, 0.5, self.wing)
        rspar = builder.spar
        pln1 = PlaneByAxes(fspar.p1, 'xz').plane
        pln2 = PlaneByAxes(fspar.p2, 'xz').plane
        builder = RibsBetweenPlanesByDistance('rib', pln1, pln2, 36., fspar,
                                              rspar, self.wing)
        self.assertEqual(builder.nribs, 12)
        self.assertAlmostEqual(builder.spacing, 35.0892, delta=0.001)
        self.assertEqual(builder.next_index, 13)
        for rib in builder.ribs:
            self.assertIsInstance(rib, Rib)

    def test_ribs_along_curve_by_number(self):
        builder = SparByParameters('fspar', 0.15, 0.15, 0.15, 0.5, self.wing)
        fspar = builder.spar
        builder = SparByParameters('rspar', 0.65, 0.15, 0.65, 0.5, self.wing)
        rspar = builder.spar
        root = RibByPoints('root', fspar.p1, rspar.p1, self.wing).rib
        tip = RibByPoints('root', fspar.p2, rspar.p2, self.wing).rib
        p1 = root.point_on_cref(0.5 * (root.cref.u1 + root.cref.u2))
        p2 = tip.point_on_cref(0.5 * (tip.cref.u1 + tip.cref.u2))
        curve = NurbsCurveByPoints([p1, p2]).curve
        builder = RibsAlongCurveByNumber('rib', curve, 10, fspar, rspar,
                                         self.wing)
        self.assertEqual(builder.nribs, 10)
        self.assertEqual(builder.next_index, 11)
        for rib in builder.ribs:
            self.assertIsInstance(rib, Rib)

    def test_ribs_along_curve_by_distance(self):
        builder = SparByParameters('fspar', 0.15, 0.15, 0.15, 0.5, self.wing)
        fspar = builder.spar
        builder = SparByParameters('rspar', 0.65, 0.15, 0.65, 0.5, self.wing)
        rspar = builder.spar
        root = RibByPoints('root', fspar.p1, rspar.p1, self.wing).rib
        tip = RibByPoints('root', fspar.p2, rspar.p2, self.wing).rib
        p1 = root.point_on_cref(0.5 * (root.cref.u1 + root.cref.u2))
        p2 = tip.point_on_cref(0.5 * (tip.cref.u1 + tip.cref.u2))
        curve = NurbsCurveByPoints([p1, p2]).curve
        builder = RibsAlongCurveByDistance('rib', curve, 36., fspar, rspar,
                                           self.wing)
        self.assertEqual(builder.nribs, 16)
        self.assertEqual(builder.next_index, 17)
        self.assertAlmostEqual(builder.spacing, 34.2162, delta=0.001)
        for rib in builder.ribs:
            self.assertIsInstance(rib, Rib)

    def test_bulkhead_by_surface(self):
        pln = PlaneByAxes((600., 0., 0.), 'yz').plane
        builder = BulkheadBySurface('bulkhead', pln, self.fuselage)
        bh = builder.bulkhead
        self.assertIsInstance(bh, Bulkhead)

    def test_floor_by_surface(self):
        pln = PlaneByAxes((0., 0., 0.), 'xy').plane
        builder = FloorBySurface('floor', pln, self.fuselage)
        floor = builder.floor
        self.assertIsInstance(floor, Floor)

    def test_frame_by_plane(self):
        pln = PlaneByAxes((600., 0., 0.), 'yz').plane
        builder = FrameByPlane('frame', pln, self.fuselage, 3.)
        frame = builder.frame
        self.assertIsInstance(frame, Frame)

    def test_frames_between_planes_by_number(self):
        pln1 = PlaneByAxes((600., 0., 0.), 'yz').plane
        pln2 = PlaneByAxes((800., 0., 0.), 'yz').plane
        builder = FramesBetweenPlanesByNumber('frame', pln1, pln2, 10,
                                              self.fuselage, 3.)
        self.assertEqual(builder.nframes, 10)
        self.assertEqual(builder.next_index, 11)
        self.assertAlmostEqual(builder.spacing, 18.1818, delta=0.001)

    def test_frames_between_planes_by_distance(self):
        pln1 = PlaneByAxes((600., 0., 0.), 'yz').plane
        pln2 = PlaneByAxes((800., 0., 0.), 'yz').plane
        builder = FramesBetweenPlanesByDistance('frame', pln1, pln2, 24.,
                                                self.fuselage, 3.)
        self.assertEqual(builder.nframes, 8)
        self.assertEqual(builder.next_index, 9)
        self.assertAlmostEqual(builder.spacing, 22.2222, delta=0.001)

    def test_frames_at_planes(self):
        pln1 = PlaneByAxes((600., 0., 0.), 'yz').plane
        pln2 = PlaneByAxes((700., 0., 0.), 'yz').plane
        pln3 = PlaneByAxes((800., 0., 0.), 'yz').plane
        builder = FramesByPlanes('frame', [pln1, pln2, pln3], self.fuselage,
                                 3.)
        self.assertEqual(builder.nframes, 3)
        self.assertEqual(builder.next_index, 4)
        self.assertAlmostEqual(builder.spacing, 100.0000, delta=0.001)

    def test_skin_by_solid(self):
        skin = SkinBySolid('skin', self.wing).skin
        self.assertIsInstance(skin, Skin)

    def test_skin_by_body(self):
        skin = SkinByBody('skin', self.fuselage).skin
        self.assertIsInstance(skin, Skin)


if __name__ == '__main__':
    unittest.main()
