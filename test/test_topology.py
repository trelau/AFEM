import unittest

from OCCT.TopoDS import TopoDS_Face, TopoDS_Shell, TopoDS_Solid

from afem.exchange import brep
from afem.geometry import Point
from afem.graphics import Viewer
from afem.topology import *


def show_shapes(*shapes):
    v = Viewer()
    v.add(*shapes)
    v.start()


class TestCreate(unittest.TestCase):
    """
    Test cases for topology creation.
    """

    def test_box_by_size(self):
        builder = BoxBySize(10., 10., 10.)
        self.assertIsInstance(builder.shell, TopoDS_Shell)
        self.assertIsInstance(builder.solid, TopoDS_Solid)
        self.assertIsInstance(builder.bottom_face, TopoDS_Face)
        self.assertIsInstance(builder.back_face, TopoDS_Face)
        self.assertIsInstance(builder.front_face, TopoDS_Face)
        self.assertIsInstance(builder.left_face, TopoDS_Face)
        self.assertIsInstance(builder.right_face, TopoDS_Face)
        self.assertIsInstance(builder.top_face, TopoDS_Face)

    def test_box_by_2_points(self):
        builder = BoxBy2Points((0, 0, 0), (10, 10, 10))
        self.assertIsInstance(builder.shell, TopoDS_Shell)
        self.assertIsInstance(builder.solid, TopoDS_Solid)
        self.assertIsInstance(builder.bottom_face, TopoDS_Face)
        self.assertIsInstance(builder.back_face, TopoDS_Face)
        self.assertIsInstance(builder.front_face, TopoDS_Face)
        self.assertIsInstance(builder.left_face, TopoDS_Face)
        self.assertIsInstance(builder.right_face, TopoDS_Face)
        self.assertIsInstance(builder.top_face, TopoDS_Face)

    def test_cylinder_by_axis(self):
        builder = CylinderByAxis(1, 10)
        self.assertIsInstance(builder.face, TopoDS_Face)
        self.assertIsInstance(builder.shell, TopoDS_Shell)
        self.assertIsInstance(builder.solid, TopoDS_Solid)

    def test_sphere_by_3_points(self):
        p1 = Point(0, 0, 0)
        p2 = Point(1, 0, 0)
        p3 = Point(0.5, 0.5, 0)
        builder = SphereBy3Points(p1, p2, p3)
        self.assertIsInstance(builder.face, TopoDS_Face)
        self.assertIsInstance(builder.shell, TopoDS_Shell)
        self.assertIsInstance(builder.solid, TopoDS_Solid)


class TestBop(unittest.TestCase):
    """
    Test cases for Boolean operations.
    """

    @unittest.expectedFailure
    def test_section_fail_08292017(self):
        """
        Intersection between a planar face and an edge. Should result in a
        vertex but finds nothing. Appears to be OCCT bug. Same behavior in
        SALOME 8.2.0.
        """
        shape1 = brep.read_brep('./test_io/section_fail_08292017_shape1.brep')
        shape2 = brep.read_brep('./test_io/section_fail_08292017_shape2.brep')

        section = IntersectShapes(shape1, shape2)
        self.assertTrue(section.is_done)

        verts = section.vertices
        self.assertEqual(len(verts), 1)

    def test_section_fail_0912017(self):
        """
        Intersection between wing solid and planar face near a wing cross
        section.
        """
        shape1 = brep.read_brep('./test_io/section_fail_09122017_shape1.brep')
        shape2 = brep.read_brep('./test_io/section_fail_09122017_shape2.brep')

        section = IntersectShapes(shape1, shape2)
        self.assertTrue(section.is_done)

        shape = section.shape
        self.assertTrue(CheckShape(shape).is_valid)

        self.assertEqual(len(section.edges), 6)

        builder = WiresByShape(shape)
        self.assertEqual(builder.nwires, 1)
        wire = builder.wires[0]
        self.assertTrue(wire.Closed())

    def test_fuse_fail_09142017(self):
        """
        In OCE 0.18 these shapes fail to fuse. In SALOME 8.2.0, they will fuse
        with the "Remove extra edges" set to OFF.
        """
        shape1 = brep.read_brep('./test_io/fuse_fail_09142017_shape1.brep')
        shape2 = brep.read_brep('./test_io/fuse_fail_09142017_shape2.brep')

        fuse = FuseShapes(shape1, shape2)
        self.assertTrue(fuse.is_done)

        shape = fuse.shape
        self.assertTrue(CheckShape(shape).is_valid)
        faces = ExploreShape.get_faces(shape)
        self.assertEqual(len(faces), 4)

        section_edges = fuse.section_edges
        self.assertEqual(len(section_edges), 1)

    def test_cut_fail_10232017(self):
        shape1 = brep.read_brep('./test_io/cut_fail_10232017_shape1.brep')
        shape2 = brep.read_brep('./test_io/cut_fail_10232017_shape2.brep')

        cut = CutShapes(shape1, shape2)
        self.assertTrue(cut.is_done)

        shape = cut.shape
        self.assertTrue(CheckShape(shape).is_valid)

    @unittest.expectedFailure
    def test_common_fail_02282018(self):
        """
        Common operation between a wing solid and a basis shell. The result
        is empty but should not be.
        """
        shape1 = brep.read_brep('./test_io/common_fail_02282018_shape1.brep')
        shape2 = brep.read_brep('./test_io/common_fail_02282018_shape2.brep')

        common = CommonShapes(shape1, shape2)
        self.assertTrue(common.is_done)

        shape = common.shape
        self.assertTrue(CheckShape(shape).is_valid)

        faces = ExploreShape.get_faces(shape)
        self.assertGreater(len(faces), 0)

    @unittest.expectedFailure
    def test_section_fail_03012018(self):
        """
        Intersection between xz-plane and a face. The intersection should be
        a single edge, but it's sometimes two edges with a gap in between or
        an edge and a single vertex.
        """
        shape1 = brep.read_brep('./test_io/section_fail_03012018_shape1.brep')
        shape2 = brep.read_brep('./test_io/section_fail_03012018_shape2.brep')

        section = IntersectShapes(shape1, shape2)
        self.assertTrue(section.is_done)

        shape = section.shape
        self.assertTrue(CheckShape(shape).is_valid)

        self.assertEqual(len(section.edges), 1)
        self.assertEqual(len(section.vertices), 2)


if __name__ == '__main__':
    unittest.main()
