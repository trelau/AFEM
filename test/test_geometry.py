# This file is part of AFEM which provides an engineering toolkit for airframe
# finite element modeling during conceptual design.
#
# Copyright (C) 2016-2018  Laughlin Research, LLC (info@laughlinresearch.com)
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
import unittest

from afem.geometry import *


class TestGeometryCreate(unittest.TestCase):
    """
    Test cases for geometry creation.
    """

    def test_point_by_xyz(self):
        p = PointByXYZ(1., 2., 3.).point
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, 2.)
        self.assertAlmostEqual(p.z, 3.)

    def test_point_by_array(self):
        p = PointByArray([1., 2., 3.]).point
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, 2.)
        self.assertAlmostEqual(p.z, 3.)

    def test_point_from_parameter(self):
        p1 = Point()
        p2 = Point(10., 0., 0.)
        line = LineByPoints(p1, p2).line
        builder = PointFromParameter(line, 0., 1.)
        p, u = builder.point, builder.parameter
        self.assertIsInstance(p, Point)
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, 0.)
        self.assertAlmostEqual(p.z, 0.)
        self.assertAlmostEqual(u, 1.)

    def test_points_along_curve_by_number(self):
        p1 = Point()
        p2 = Point(10., 0., 0.)
        line = LineByPoints(p1, p2).line
        builder = PointsAlongCurveByNumber(line, 3, 0., 10.)
        self.assertEqual(builder.npts, 3)
        self.assertAlmostEqual(builder.spacing, 5.)
        p1, p2, p3 = builder.points
        u1, u2, u3 = builder.parameters
        self.assertAlmostEqual(p1.x, 0.)
        self.assertAlmostEqual(p2.x, 5.)
        self.assertAlmostEqual(p3.x, 10.)
        self.assertAlmostEqual(u1, 0.)
        self.assertAlmostEqual(u2, 5.)
        self.assertAlmostEqual(u3, 10.)

    def test_points_along_curve_by_distance(self):
        p1 = Point()
        p2 = Point(10., 0., 0.)
        line = LineByPoints(p1, p2).line
        builder = PointsAlongCurveByDistance(line, 5., 0., 10.)
        self.assertEqual(builder.npts, 3)
        self.assertAlmostEqual(builder.spacing, 5.)
        p1, p2, p3 = builder.points
        u1, u2, u3 = builder.parameters
        self.assertAlmostEqual(p1.x, 0.)
        self.assertAlmostEqual(p2.x, 5.)
        self.assertAlmostEqual(p3.x, 10.)
        self.assertAlmostEqual(u1, 0.)
        self.assertAlmostEqual(u2, 5.)
        self.assertAlmostEqual(u3, 10.)

    def test_direction_by_xyz(self):
        d = DirectionByXYZ(1., 0., 0.).direction
        self.assertAlmostEqual(d.i, 1.)
        self.assertAlmostEqual(d.j, 0.)
        self.assertAlmostEqual(d.k, 0.)

    def test_direction_by_array(self):
        d = DirectionByArray([1., 0., 0.]).direction
        self.assertAlmostEqual(d.i, 1.)
        self.assertAlmostEqual(d.j, 0.)
        self.assertAlmostEqual(d.k, 0.)

    def test_direction_by_points(self):
        p1 = Point()
        d = DirectionByPoints(p1, (10., 0., 0.)).direction
        self.assertAlmostEqual(d.i, 1.)
        self.assertAlmostEqual(d.j, 0.)
        self.assertAlmostEqual(d.k, 0.)

    def test_vector_by_xyz(self):
        v = VectorByXYZ(1., 2., 3.).vector
        self.assertAlmostEqual(v.x, 1.)
        self.assertAlmostEqual(v.y, 2.)
        self.assertAlmostEqual(v.z, 3.)

    def test_vector_by_array(self):
        v = VectorByArray([1., 2., 3.]).vector
        self.assertAlmostEqual(v.x, 1.)
        self.assertAlmostEqual(v.y, 2.)
        self.assertAlmostEqual(v.z, 3.)

    def test_vector_by_points(self):
        p1 = Point()
        v = VectorByPoints(p1, (1., 2., 3.)).vector
        self.assertAlmostEqual(v.x, 1.)
        self.assertAlmostEqual(v.y, 2.)
        self.assertAlmostEqual(v.z, 3.)

    def test_line_by_vector(self):
        line = LineByVector((0, 0, 0,), (1, 0, 0)).line
        self.assertIsInstance(line, Line)
        p = line.eval(10.)
        self.assertAlmostEqual(p.x, 10.)

    def test_line_by_points(self):
        line = LineByPoints((0, 0, 0,), (10, 0, 0)).line
        self.assertIsInstance(line, Line)
        p = line.eval(10.)
        self.assertAlmostEqual(p.x, 10.)

    def circle_by_normal(self):
        circle = CircleByNormal((0, 0, 0), (1, 0, 0), 1).circle
        self.assertIsInstance(circle, Circle)
        p = circle.eval(0)
        self.assertAlmostEqual(p.x, 0.)
        self.assertAlmostEqual(p.y, 0.)
        self.assertAlmostEqual(p.z, 1.)

    def circle_by_plane(self):
        pln = PlaneByAxes(axes='xy').plane
        circle = CircleByPlane((0, 0, 0), pln, 1).circle
        self.assertIsInstance(circle, Circle)
        p = circle.eval(0)
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, 0.)
        self.assertAlmostEqual(p.z, 0.)

    def test_circle_by_3_points(self):
        p1 = Point(0, 0, 0)
        p2 = Point(1, 0, 0)
        p3 = Point(0.5, 0.5, 0)
        circle = CircleBy3Points(p1, p2, p3).circle
        self.assertAlmostEqual(circle.radius, 0.5, 6)
        self.assertAlmostEqual(circle.center.x, 0.5)
        self.assertAlmostEqual(circle.center.y, 0.)
        self.assertAlmostEqual(circle.center.y, 0.)

    def test_nurbs_curve_by_data(self):
        cp = [(0, 0, 0), (10, 0, 0)]
        uk = [0, 1]
        m = [2, 2]
        c = NurbsCurveByData(cp, uk, m, 1).curve
        p = c.eval(0.5)
        self.assertIsInstance(c, NurbsCurve)
        self.assertAlmostEqual(p.x, 5.)

    def test_nurbs_curve_by_interp(self):
        qp = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]
        c = NurbsCurveByInterp(qp).curve
        p = c.eval(5.)
        self.assertIsInstance(c, NurbsCurve)
        self.assertEqual(c.p, 2)
        self.assertAlmostEqual(p.x, 3.536, places=3)
        self.assertAlmostEqual(p.y, 4.571, places=3)
        self.assertAlmostEqual(p.z, 0.)

    def test_nurbs_curve_by_approx(self):
        qp = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]
        c = NurbsCurveByApprox(qp).curve
        p = c.eval(0.5)
        self.assertIsInstance(c, NurbsCurve)
        self.assertEqual(c.p, 3)
        self.assertAlmostEqual(p.x, 5.)
        self.assertAlmostEqual(p.y, 5.)
        self.assertAlmostEqual(p.z, 0.)

    def test_nurbs_curve_by_points(self):
        qp = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]
        c = NurbsCurveByPoints(qp).curve
        p = c.eval(0.5)
        self.assertIsInstance(c, NurbsCurve)
        self.assertEqual(c.p, 1)
        self.assertAlmostEqual(p.x, 5.)
        self.assertAlmostEqual(p.y, 5.)
        self.assertAlmostEqual(p.z, 0.)

    def test_trimmed_curve_by_parameters(self):
        qp = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]
        basis_curve = NurbsCurveByInterp(qp).curve
        c = TrimmedCurveByParameters(basis_curve, 1., 9.).curve
        self.assertIsInstance(c, TrimmedCurve)
        self.assertAlmostEqual(c.u1, 1.)
        self.assertAlmostEqual(c.u2, 9.)

    def test_trimmed_curve_by_curve(self):
        qp = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]
        basis_curve = NurbsCurveByInterp(qp).curve
        c = TrimmedCurveByCurve(basis_curve).curve
        self.assertIsInstance(c, TrimmedCurve)
        self.assertAlmostEqual(c.u1, 0.)
        self.assertAlmostEqual(c.u2, 14.142136, places=6)

    def test_trimmed_curve_by_points(self):
        qp = [(0, 0, 0), (5, 5, 0), (10, 0, 0)]
        basis_curve = NurbsCurveByInterp(qp).curve
        p1 = basis_curve.eval(1.)
        p2 = basis_curve.eval(9.)
        c = TrimmedCurveByPoints(basis_curve, p1, p2).curve
        self.assertAlmostEqual(c.u1, 1.)
        self.assertAlmostEqual(c.u2, 9.)

    def test_plane_by_normal(self):
        pln = PlaneByNormal((0, 0, 0), (0, 0, 1)).plane
        p = pln.eval(1., 1.)
        self.assertIsInstance(pln, Plane)
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, 1.)
        self.assertAlmostEqual(p.z, 0.)

    def test_plane_by_axes(self):
        xy_pln = PlaneByAxes(axes='xy').plane
        xz_pln = PlaneByAxes(axes='xz').plane
        yz_pln = PlaneByAxes(axes='yz').plane
        self.assertIsInstance(xy_pln, Plane)
        self.assertIsInstance(xz_pln, Plane)
        self.assertIsInstance(yz_pln, Plane)
        p1 = xy_pln.eval(1., 1.)
        p2 = xz_pln.eval(1., 1.)
        p3 = yz_pln.eval(1., 1.)
        self.assertAlmostEqual(p1.x, 1.)
        self.assertAlmostEqual(p1.y, 1.)
        self.assertAlmostEqual(p1.z, 0.)
        self.assertAlmostEqual(p2.x, 1.)
        self.assertAlmostEqual(p2.y, 0.)
        self.assertAlmostEqual(p2.z, 1.)
        self.assertAlmostEqual(p3.x, 0.)
        self.assertAlmostEqual(p3.y, 1.)
        self.assertAlmostEqual(p3.z, 1.)

    def test_plane_by_points(self):
        pln = PlaneByPoints((0, 0, 0), (1, 0, 0), (0, 1, 0)).plane
        self.assertIsInstance(pln, Plane)
        p = pln.eval(1., 1.)
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, 1.)
        self.assertAlmostEqual(p.z, 0.)

    def test_plane_by_Approx(self):
        pln = PlaneByApprox([(0, 0, 0), (1, 0, 0), (0, 1, 0)]).plane
        self.assertIsInstance(pln, Plane)
        p = pln.eval(1., 1.)
        self.assertAlmostEqual(p.x, 1.333, places=3)
        self.assertAlmostEqual(p.y, 1.333, places=3)
        self.assertAlmostEqual(p.z, 0.)

    def test_plane_by_parameter(self):
        line = LineByPoints((0, 0, 0), (10, 0, 0)).line
        pln = PlaneFromParameter(line, 0., 1.).plane
        self.assertIsInstance(pln, Plane)
        p = pln.eval(1., 1.)
        self.assertAlmostEqual(p.x, 1.)
        self.assertAlmostEqual(p.y, -1.)
        self.assertAlmostEqual(p.z, 1.)


if __name__ == '__main__':
    unittest.main()
