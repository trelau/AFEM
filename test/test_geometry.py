import unittest

from afem.geometry import *


class TestCreate(unittest.TestCase):
    """
    Test cases for geometry creation.
    """

    def test_circle_by_3_points(self):
        p1 = Point(0, 0, 0)
        p2 = Point(1, 0, 0)
        p3 = Point(0.5, 0.5, 0)
        circle = CircleBy3Points(p1, p2, p3).circle
        self.assertAlmostEqual(circle.radius, 0.5, 6)
        self.assertAlmostEqual(circle.center.x, 0.5)
        self.assertAlmostEqual(circle.center.y, 0.)
        self.assertAlmostEqual(circle.center.y, 0.)


if __name__ == '__main__':
    unittest.main()
