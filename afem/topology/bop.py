#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2017 Laughlin Research, L.L.C.
#
# This file is subject to the license agreement that was delivered
# with this source code.
#
# THE SOFTWARE AND INFORMATION ARE PROVIDED ON AN "AS IS" BASIS,
# WITHOUT ANY WARRANTIES OR REPRESENTATIONS EXPRESS, IMPLIED OR
# STATUTORY; INCLUDING, WITHOUT LIMITATION, WARRANTIES OF QUALITY,
# PERFORMANCE, MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

from datetime import datetime
from warnings import warn

from OCCT.BOPAlgo import BOPAlgo_MakerVolume, BOPAlgo_Options
from OCCT.BRepAlgoAPI import (BRepAlgoAPI_Common, BRepAlgoAPI_Cut,
                              BRepAlgoAPI_Fuse, BRepAlgoAPI_Section,
                              BRepAlgoAPI_Splitter)
from OCCT.BRepFeat import BRepFeat_MakeCylindricalHole, BRepFeat_SplitShape
from OCCT.Message import Message_Gravity
from OCCT.TopTools import TopTools_SequenceOfShape
from OCCT.TopoDS import TopoDS_Face, TopoDS, TopoDS_Wire

from afem.geometry.check import CheckGeom
from afem.occ.utils import (to_lst_from_toptools_listofshape,
                            to_toptools_listofshape)
from afem.topology.check import CheckShape
from afem.topology.explore import ExploreShape, ExploreWire

__all__ = ["BopCore", "BopAlgo", "FuseShapes", "CutShapes", "CommonShapes",
           "IntersectShapes", "SplitShapes", "VolumesFromShapes",
           "CutCylindricalHole", "LocalSplit", "SplitShapeByEdges",
           "TrimOpenWire"]

# Turn on parallel Boolean execution by default
BOPAlgo_Options.SetParallelMode_(True)

# Message gravities
_gravities = [Message_Gravity.Message_Trace, Message_Gravity.Message_Info,
              Message_Gravity.Message_Warning, Message_Gravity.Message_Alarm,
              Message_Gravity.Message_Fail]


class BopCore(object):
    """
    Core class for Boolean operations and enabling attributes and methods for
    rebuilding shapes.
    """

    def __init__(self):
        self._bop = None

    def build(self):
        """
        Build the results.

        :return: None.
        """
        if isinstance(self._bop, BOPAlgo_MakerVolume):
            self._bop.Perform()
        else:
            self._bop.Build()

    @property
    def is_done(self):
        """
        :return: *True* if operation is done, *False* if not.
        :rtype: bool
        """
        if isinstance(self._bop, (BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            return not self._bop.HasErrors()
        return self._bop.IsDone()

    @property
    def shape(self):
        """
        :return: The resulting shape.
        :rtype: OCCT.TopoDS.TopoDS_Shape
        """
        return self._bop.Shape()

    def modified(self, shape):
        """
        Return a list of shapes modified from the given shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: List of modified shapes.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._bop.Modified(shape))

    def generated(self, shape):
        """
        Return a list of shapes generated from the given shape.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: List of generated shapes.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._bop.Generated(shape))

    def is_deleted(self, shape):
        """
        Check to see if shape is deleted.

        :param OCCT.TopoDS.TopoDS_Shape shape: The shape.

        :return: *True* if deleted, *False* if not.
        :rtype: bool
        """
        return self._bop.IsDeleted(shape)


class BopAlgo(BopCore):
    """
    Base class for Boolean operations.

    :param shape1: The first shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or None
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.
    :param bop: The OpenCASCADE class for the Boolean operation.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.
    """

    def __init__(self, shape1, shape2, fuzzy_val, nondestructive, bop):
        super(BopAlgo, self).__init__()

        self._bop = bop()

        if fuzzy_val is not None:
            self._bop.SetFuzzyValue(fuzzy_val)

        if nondestructive:
            self._bop.SetNonDestructive(True)
        else:
            self._bop.SetNonDestructive(False)

        if CheckShape.is_shape(shape1) and CheckShape.is_shape(shape2):
            self.set_args([shape1])
            self.set_tools([shape2])
            self.build()

    @staticmethod
    def set_parallel_mode(flag):
        """
        Global option to set the Boolean operations for parallel execution.

        :param bool flag: Option for parallel execution. *True* turns parallel
            execution on, *False* turns it off.

        :return: None.
        """
        BOPAlgo_Options.SetParallelMode_(flag)

    def debug(self, path='.'):
        """
        Export files for debugging Boolean operations.

        :param path:

        :return:
        """
        # Generate a suffix using timestamp
        now = datetime.now()
        timestamp = str(now.timestamp())

        # Operation name for prefix
        op = str(self.__class__.__name__)

        # Info file
        fn = ''.join([path, '/', op, '.info.', timestamp, '.txt'])
        info = open(fn, 'w')
        info.write('Date (M-D-Y): {}-{}-{}\n'.format(now.month, now.day,
                                                     now.year))
        info.write('Operation: {}\n'.format(op))
        info.write('Parallel: {}\n'.format(self._bop.RunParallel()))
        info.write('Fuzzy value: {}\n'.format(self._bop.FuzzyValue()))
        info.write('Nondestructive: {}\n'.format(self._bop.NonDestructive()))

        # Errors and warnings report
        msg_report = self._bop.GetReport()
        for gravity in _gravities:
            msg_list = msg_report.GetAlerts(gravity)
            if msg_list.Size() == 0:
                continue
            info.write('Messages:\n')
            for msg in msg_list:
                info.write('\t{}\n'.format(msg.GetMessageKey()))

        # Avoid circular imports
        from afem.exchange.brep import write_brep
        from afem.topology.create import CompoundByShapes

        # Arguments
        args = self.arguments
        if args:
            shape1 = CompoundByShapes(args).compound
            fn = ''.join([path, '/', op, '.shape1.', timestamp, '.brep'])
            write_brep(shape1, fn)

        # Tools
        tools = self.tools
        if tools:
            shape2 = CompoundByShapes(tools).compound
            fn = ''.join([path, '/', op, '.shape2.', timestamp, '.brep'])
            write_brep(shape2, fn)

    @property
    def arguments(self):
        """
        :return: The arguments.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._bop.Arguments())

    @property
    def tools(self):
        """
        :return: The tools.
        :rtype: list[OCCT.TopoDS.TopoDS_Shape]
        """
        return to_lst_from_toptools_listofshape(self._bop.Tools())

    def set_args(self, shapes):
        """
        Set the arguments.

        :param list[OCCT.TopoDS.TopoDS_Shape] shapes: The arguments.

        :return: None.
        """
        if isinstance(self._bop, BOPAlgo_MakerVolume):
            for shape in shapes:
                self._bop.AddArgument(shape)
            return None
        args = to_toptools_listofshape(shapes)
        self._bop.SetArguments(args)

    def set_tools(self, shapes):
        """
        Set the tools.

        :param list[OCCT.TopoDS.TopoDS_Shape] shapes: The tools.

        :return: None.
        """
        if isinstance(self._bop, BOPAlgo_MakerVolume):
            warn('Setting tools not available. Doing nothing.', RuntimeWarning)
            return None

        tools = to_toptools_listofshape(shapes)
        self._bop.SetTools(tools)

    @property
    def vertices(self):
        """
        :return: The vertices of the resulting shape.
        :rtype: list[OCCT.TopoDS.TopoDS_Vertex]
        """
        return ExploreShape.get_vertices(self.shape)

    @property
    def edges(self):
        """
        :return: The edges of the resulting shape.
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        return ExploreShape.get_edges(self.shape)

    def refine_edges(self):
        """
        Fuse C1 edges.

        :return: None.
        """
        if isinstance(self._bop, (BRepAlgoAPI_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            warn('Refining edges not available. Doing nothing.',
                 RuntimeWarning)
        else:
            self._bop.RefineEdges()

    @property
    def fuse_edges(self):
        """
        :return: The result flag of edge refining.
        :rtype: bool
        """
        if isinstance(self._bop, (BRepAlgoAPI_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            return False
        else:
            return self._bop.FuseEdges()

    @property
    def section_edges(self):
        """
        :return: A list of section edges as a result of intersection between
            the shapes.
        :rtype: list[OCCT.TopoDS.TopoDS_Edge]
        """
        if isinstance(self._bop, (BRepAlgoAPI_Splitter, BOPAlgo_MakerVolume,
                                  BRepFeat_MakeCylindricalHole)):
            warn('Getting section edges not available. Returning an empty '
                 'list.', RuntimeWarning)
            return []
        else:
            return to_lst_from_toptools_listofshape(self._bop.SectionEdges())

    @property
    def has_modified(self):
        """
        :return: *True* if there is at least one modified shape.
        :rtype: bool
        """
        return self._bop.HasModified()

    @property
    def has_generated(self):
        """
        :return: *True* if there is at least one generated shape.
        :rtype: bool
        """
        return self._bop.HasGenerated()

    @property
    def has_deleted(self):
        """
        :return: *True* if there is at least one deleted shape.
        :rtype: bool
        """
        return self._bop.HasDeleted()


class FuseShapes(BopAlgo):
    """
    Boolean fuse operation.

    :param shape1: The first shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or None
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Fuse_.

    .. _BRepAlgoAPI_Fuse: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_b_rep_algo_a_p_i___fuse.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = FuseShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = FuseShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, fuzzy_val=None,
                 nondestructive=False):
        super(FuseShapes, self).__init__(shape1, shape2, fuzzy_val,
                                         nondestructive, BRepAlgoAPI_Fuse)


class CutShapes(BopAlgo):
    """
    Boolean cut operation.

    :param shape1: The first shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or None
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Cut_.

    .. _BRepAlgoAPI_Cut: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_b_rep_algo_a_p_i___cut.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 0., 0.), (6., 0., 0.)).edge
    >>> bop = CutShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = CutShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, fuzzy_val=None,
                 nondestructive=False):
        super(CutShapes, self).__init__(shape1, shape2, fuzzy_val,
                                        nondestructive, BRepAlgoAPI_Cut)


class CommonShapes(BopAlgo):
    """
    Boolean common operation.

    :param shape1: The first shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or None
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Common_.

    .. _BRepAlgoAPI_Common: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_b_rep_algo_a_p_i___common.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 0., 0.), (6., 0., 0.)).edge
    >>> bop = CommonShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = CommonShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, fuzzy_val=None,
                 nondestructive=False):
        super(CommonShapes, self).__init__(shape1, shape2, fuzzy_val,
                                           nondestructive, BRepAlgoAPI_Common)


class IntersectShapes(BopAlgo):
    """
    Boolean intersect operation.

    :param shape1: The first shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param shape2: The second shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param bool compute_pcurve1: Option to compute p-curves on shape 1.
    :param bool compute_pcurve2: Option to compute p-curves on shape 2.
    :param bool approximate: Option to approximate intersection curves.
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    For more information see BRepAlgoAPI_Section_.

    .. _BRepAlgoAPI_Section: https://www.opencascade.com/doc/occt-7.2.0/refman/html/class_b_rep_algo_a_p_i___section.html

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = IntersectShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = IntersectShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, compute_pcurve1=False,
                 compute_pcurve2=False, approximate=False, fuzzy_val=None,
                 nondestructive=False):
        super(IntersectShapes, self).__init__(None, None, fuzzy_val,
                                              nondestructive,
                                              BRepAlgoAPI_Section)

        self._bop.ComputePCurveOn1(compute_pcurve1)
        self._bop.ComputePCurveOn2(compute_pcurve2)
        self._bop.Approximation(approximate)

        build1, build2 = False, False
        if CheckShape.is_shape(shape1):
            self._bop.Init1(shape1)
            build1 = True
        elif CheckGeom.is_surface(shape1):
            self._bop.Init1(shape1.object)
            build1 = True

        if CheckShape.is_shape(shape2):
            self._bop.Init2(shape2)
            build2 = True
        elif CheckGeom.is_surface(shape2):
            self._bop.Init2(shape2.object)
            build2 = True

        if build1 and build2:
            self._bop.Build()

    def has_ancestor_face1(self, edge):
        """
        Get the ancestor face on the intersection edge on the first shape
        if available.

        :param OCCT.TopoDS.TopoDS_Edge edge: The edge.

        :return: *True* and the face if available, *False* and *None* if not.
        :rtype: tuple(bool, OCCT.TopoDS.TopoDS_Face or None)
        """
        f = TopoDS_Face()
        if self._bop.HasAncestorFaceOn1(edge, f):
            return True, f
        return False, None

    def has_ancestor_face2(self, edge):
        """
        Get the ancestor face on the intersection edge on the second shape
        if available.

        :param OCCT.TopoDS.TopoDS_Edge edge: The edge.

        :return: *True* and the face if available, *False* and *None* if not.
        :rtype: tuple(bool, OCCT.TopoDS.TopoDS_Face or None)
        """
        f = TopoDS_Face()
        if self._bop.HasAncestorFaceOn2(edge, f):
            return True, f
        return False, None


class SplitShapes(BopAlgo):
    """
    Split arbitrary shapes. This is a wrapper for the SALOME
    GEOMAlgo_Splitter tool.

    :param shape1: The first shape.
    :type shape1: OCCT.TopoDS.TopoDS_Shape or None
    :param shape2: The second shape.
    :type shape2: OCCT.TopoDS.TopoDS_Shape or None
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    .. note::

        If *shape1* or *shape2* is *None* then the user is expected to manually
        set the arguments and tools and build the result.

    Usage:

    >>> from afem.topology import *
    >>> e1 = EdgeByPoints((0., 0., 0.), (10., 0., 0.)).edge
    >>> e2 = EdgeByPoints((5., 1., 0.), (5., -1., 0.)).edge
    >>> bop = SplitShapes(e1, e2)
    >>> assert bop.is_done
    >>> shape = bop.shape
    >>> # Setting arguments and tools
    >>> bop = SplitShapes()
    >>> bop.set_args([e1])
    >>> bop.set_tools([e2])
    >>> bop.build()
    >>> assert bop.is_done
    """

    def __init__(self, shape1=None, shape2=None, fuzzy_val=None,
                 nondestructive=False):
        super(SplitShapes, self).__init__(shape1, shape2, fuzzy_val,
                                          nondestructive, BRepAlgoAPI_Splitter)


class VolumesFromShapes(BopAlgo):
    """
    Build solids from a list of shapes.

    :param list[OCCT.TopoDS.TopoDS_Shape] shapes: The shapes.
    :param bool intersect: Option to intersect the shapes before building
        solids.
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.
    """

    def __init__(self, shapes, intersect=False, fuzzy_val=None,
                 nondestructive=False):
        super(VolumesFromShapes, self).__init__(None, None, fuzzy_val,
                                                nondestructive,
                                                BOPAlgo_MakerVolume)

        self.set_args(shapes)

        if intersect:
            self._bop.SetIntersect(True)
        else:
            self._bop.SetIntersect(False)

        self.build()

        self._solids = []
        for solid in ExploreShape.get_solids(self.shape):
            self._solids.append(TopoDS.Solid_(solid))

    @property
    def box(self):
        """
        :return: The bounding box of all provided shapes.
        :rtype: OCCT.TopoDS.TopoDS_Solid
        """
        return self._bop.Box()

    @property
    def nsolids(self):
        """
        :return: The number of solids in the shape.
        :rtype: int
        """
        return len(self.solids)

    @property
    def solids(self):
        """
        :return: The list of solids.
        :rtype: list[OCCT.TopoDS.TopoDS_Solid]
        """
        return self._solids


class CutCylindricalHole(BopAlgo):
    """
    Cut a cylindrical hole on a shape.

    :param shape: The shape.
    :type shape: OCCT.TopoDS.TopoDS_Shape
    :param float radius: The radius of the hole.
    :param afem.geometry.entities.Axis1: The axis for the hole.
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln = PlaneByAxes().plane
    >>> face = FaceByPlane(pln, -2., 2., -2., 2.).face
    >>> # The small offset is a workaround for planar cuts
    >>> p = Point(0., 0.1, 0.)
    >>> d = Direction(0., 1., 0.)
    >>> ax1 = Axis1(p, d)
    >>> bop = CutCylindricalHole(face, 1., ax1)
    >>> assert bop.is_done
    """

    def __init__(self, shape, radius, ax1, fuzzy_val=None,
                 nondestructive=False):
        super(CutCylindricalHole, self).__init__(None, None, fuzzy_val,
                                                 nondestructive,
                                                 BRepFeat_MakeCylindricalHole)

        self._bop.Init(shape, ax1)
        self._bop.Perform(radius)


class LocalSplit(BopCore):
    """
    Perform a local split of a shape in the context of a basis shape. This tool
    only splits faces.

    :param OCCT.TopoDS.TopoDS_Shape shape: The local shape.
    :param tool: The tool to split with.
    :type tool: OCCT.TopoDS.TopoDS_Shape or afem.geometry.entities.Surface
    :param OCCT.TopoDS.TopoDS_Shape basis_shape: The basis shape that the local
        shape is part of.
    :param bool approximate: Option to approximate intersection curves.
    :param float fuzzy_val: Fuzzy tolerance value.
    :param bool nondestructive: Option to not modify the input shapes.

    Usage:

    >>> from afem.geometry import *
    >>> from afem.topology import *
    >>> pln = PlaneByAxes().plane
    >>> builder = SolidByPlane(pln, 5., 5., 5.)
    >>> box = builder.solid
    >>> tool = PlaneByAxes(axes='xy').plane
    >>> face = ExploreShape.get_faces(box)[0]
    >>> split = LocalSplit(face, tool, box)
    >>> assert split.is_done
    """

    def __init__(self, shape, tool, basis_shape, approximate=False,
                 fuzzy_val=None, nondestructive=False):
        super(LocalSplit, self).__init__()

        # Intersect
        section = IntersectShapes(shape, tool, True, False, approximate,
                                  fuzzy_val, nondestructive)
        sec_edges = section.edges

        # Split
        self._bop = BRepFeat_SplitShape(basis_shape)
        for e in sec_edges:
            status, f = section.has_ancestor_face1(e)
            if status:
                self._bop.Add(e, f)
        self.build()


class SplitShapeByEdges(BopCore):
    """
    Split a shape using edges.

    :param OCCT.TopoDS.TopoDS_Shape shape: The basis shape.
    :param edges: The edges to split the shape with. If provided, then the
        results will be built during initialization. If none are provided then
        the user is expected to add edges and build manually.
    :type edges: collections.Sequence(OCCT.TopoDS.TopoDS_Edge) or None
    :param bool check_interior: Option to check internal intersections.

    """

    def __init__(self, shape, edges=None, check_interior=True):
        super(SplitShapeByEdges, self).__init__()

        self._bop = BRepFeat_SplitShape(shape)

        if not check_interior:
            self._bop.SetCheckInterior(False)

        if edges is not None:
            edge_seq = TopTools_SequenceOfShape()
            for e in edges:
                edge_seq.Append(e)
            self._bop.Add(edge_seq)
            self.build()

    def add_edges(self, shapes):
        """
        Add splittings edges or wires for the initial shape.

        :param collections.Sequence(OCCT.TopoDS.TopoDS_Shape) shapes: The
            splitting edges or wires.

        :return: *True* if added, *False* if not.
        :rtype: bool
        """
        seq = TopTools_SequenceOfShape()
        for s in shapes:
            seq.Append(s)
        return self._bop.Add(seq)

    def add_wire_on_face(self, w, f):
        """
        Add the wire on the face.

        :param OCCT.TopoDS.TopoDS_Wire w: The wire.
        :param OCCT.TopoDS.TopoDS_Face f: The face.

        :return: None.
        """
        self._bop.Add(w, f)

    def add_edge_on_face(self, e, f):
        """
        Add the edge on the face.

        :param OCCT.TopoDS.TopoDS_Edge e: The edge.
        :param OCCT.TopoDS.TopoDS_Face f: The face.

        :return: None
        """
        self._bop.Add(e, f)

    def add_edges_on_face(self, e, f):
        """

        :param collections.Sequence(OCCT.TopoDS.TopoDS_Edge) e: The edges.
        :param OCCT.TopoDS.TopoDS_Face f: The face.

        :return: None
        """
        # Avoid circular imports
        from afem.topology.create import CompoundByShapes

        cmp = CompoundByShapes(e).compound
        self._bop.Add(cmp, f)

    def add_edge_on_edge(self, e1, e2):
        """
        Add the edge on an existing edge.

        :param OCCT.TopoDS.TopoDS_Edge e1: The edge.
        :param OCCT.TopoDS.TopoDS_Edge e2: The existing edge.

        :return: None.
        """
        self._bop.Add(e1, e2)


class TrimOpenWire(object):
    """
    Trim an open wire between two shapes.
    """

    def __init__(self, wire, shape1, shape2):
        wire = CheckShape.to_wire(wire)
        if not isinstance(wire, TopoDS_Wire):
            msg = 'Invalid type for wire.'
            raise TypeError(msg)

        if wire.Closed():
            msg = 'Closed wires are not supported.'
            raise TypeError(msg)

        shape1 = CheckShape.to_shape(shape1)
        shape2 = CheckShape.to_shape(shape2)

        # Split wire with shapes
        from afem.topology.create import CompoundByShapes

        other_shape = CompoundByShapes([shape1, shape2]).compound
        split = SplitShapes(wire, other_shape)
        split_wire = ExploreShape.get_wires(split.shape)[0]

        # Get new vertices
        old_verts = ExploreWire(wire).ordered_vertices
        wire_exp = ExploreWire(split_wire)
        all_verts = wire_exp.ordered_vertices
        new_verts = [v for v in all_verts if v not in old_verts]
        if len(new_verts) > 2:
            msg = 'More than two split locations is not supported.'
            raise RuntimeError(msg)

        # Find index of new vertices and use that to extract edges
        i1 = all_verts.index(new_verts[0])
        i2 = all_verts.index(new_verts[1])
        ordered_edges = wire_exp.edges
        first_edges = ordered_edges[:i1]
        trimmed_edges = ordered_edges[i1:i2]
        last_edges = ordered_edges[i2:]

        # Avoid circular imports
        from afem.topology.create import WireByEdges

        # Collect data and build trimmed wires
        self._first_wire = None
        self._trimmed_wire = None
        self._last_wire = None

        self._split_wire = split_wire
        if len(first_edges) > 0:
            self._first_wire = WireByEdges(*first_edges).wire
        if len(trimmed_edges) > 0:
            self._trimmed_wire = WireByEdges(*trimmed_edges).wire
        if len(last_edges) > 0:
            self._last_wire = WireByEdges(*last_edges).wire
        self._new_verts = new_verts

    @property
    def split_wire(self):
        """
        :return: The wire after splitting.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._split_wire

    @property
    def first_wire(self):
        """
        :return: The first trimmed segment.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._first_wire

    @property
    def last_wire(self):
        """
        :return: The last trimmed segment.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._last_wire

    @property
    def trimmed_wire(self):
        """
        :return: The interior trimmed segment.
        :rtype: OCCT.TopoDS.TopoDS_Wire
        """
        return self._trimmed_wire

    @property
    def new_vertices(self):
        """
        :return: New vertices after splitting.
        :rtype: list(OCCT.TopoDS.TopoDS_Vertex)
        """
        return self._new_verts


if __name__ == "__main__":
    import doctest

    doctest.testmod()
