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

import json

from OCCT.BRep import BRep_Builder, BRep_Tool
from OCCT.BRepBuilderAPI import (BRepBuilderAPI_MakeFace,
                                 BRepBuilderAPI_MakeWire, BRepBuilderAPI_Sewing)
from OCCT.BRepCheck import BRepCheck_Analyzer, BRepCheck_NoError
from OCCT.BRepGProp import BRepGProp
from OCCT.GCPnts import GCPnts_QuasiUniformDeflection
from OCCT.GProp import GProp_GProps
from OCCT.Geom import Geom_Plane
from OCCT.GeomAdaptor import GeomAdaptor_Curve
from OCCT.GeomLib import GeomLib_IsPlanarSurface
from OCCT.IFSelect import (IFSelect_ItemsByEntity, IFSelect_RetDone,
                           IFSelect_RetVoid)
from OCCT.Interface import Interface_Static
from OCCT.STEPControl import STEPControl_Reader
from OCCT.ShapeAnalysis import ShapeAnalysis
from OCCT.ShapeFix import ShapeFix_Solid, ShapeFix_Wire
from OCCT.ShapeUpgrade import (ShapeUpgrade_ShapeDivideClosed,
                               ShapeUpgrade_SplitSurface,
                               ShapeUpgrade_UnifySameDomain)
from OCCT.TColStd import TColStd_HSequenceOfReal
from OCCT.TopAbs import TopAbs_COMPOUND, TopAbs_FACE
from OCCT.TopExp import TopExp_Explorer
from OCCT.TopoDS import TopoDS_Compound, TopoDS_Iterator, TopoDS_Shell

from afem.config import Settings, logger
from afem.geometry.create import (PointFromParameter, NurbsSurfaceByInterp,
                                  NurbsCurveByPoints, NurbsCurveByApprox)
from afem.geometry.entities import NurbsCurve, NurbsSurface
from afem.geometry.utils import chord_parameters
from afem.oml.entities import Body
from afem.topology.check import CheckShape
from afem.topology.create import (CompoundByShapes, FaceBySurface, EdgeByCurve,
                                  WireByEdges)
from afem.topology.explore import ExploreShape
from afem.topology.fix import FixShape
from afem.topology.modify import ShapeBSplineRestriction
from afem.topology.props import LinearProps

__all__ = ["ImportVSP"]


class ImportVSP(object):
    """
    Tool for importing, pre-processing, and translating OpenVSP models.

    :param str fn: The full path to a file to import. If provided, the file
        extension will be used to guess which import method to call. For
        example, if the file ended with 'stp' or 'step', then the
        *import_step()* method is called during initialization.
    :param bool divide_closed: Option to divide closed faces.
    :param bool bspline_restrict: Option to attempt to refit the underlying
        OpenVSP surfaces with surfaces of higher continuity. This may take a
        long time to run and no guarantee it will be successful. At the moment
        this is only used on known wing components. This method is experimental.
    :param bool reloft: For wings that are not split, this option will
        extract isocurves at each spanwise cross section, tessellate the
        curve using *tol*, and then approximate the section using a C1
        continuous curve. These curves are then used to generate a new wing
        surface. This method is experimental.
    :param float tol: Tolerance for approximation if *bspline_restrict* or
        *reloft* is *True*.

    :raise TypeError: If a file is provided but the extension is not recognized
        or supported.
    """

    def __init__(self, fn=None, divide_closed=True, bspline_restrict=False,
                 reloft=False, tol=0.01):
        self._bodies = {}
        self._divide = divide_closed
        self._restrict = bspline_restrict
        self._reloft = reloft
        self._tol = tol

        if fn is not None:
            if fn.endswith('.step') or fn.endswith('.stp'):
                self.import_step(fn)
            else:
                raise TypeError('File extension not supported.')

    @property
    def bodies(self):
        """
        :return: List of Body instances.
        :rtype: list[afem.oml.entities.Body]
        """
        return list(self._bodies.values())

    def clear(self):
        """
        Clear any translated data.

        :return: None.
        """
        self._bodies.clear()

    def get_body(self, name):
        """
        Return Body by name.

        :param str name: Body name.

        :return: OML body.
        :rtype: afem.oml.entities.Body

        :raise KeyError: If *name* is not present.
        """
        return self._bodies[name]

    def get_bodies(self):
        """
        Return all Body instances in a list.

        :return: List of Body instances.
        :rtype: list[afem.oml.entities.Body]
        """
        return list(self._bodies.values())

    def import_step(self, fn):
        """
        Import a STEP file generated by the OpenVSP version that has been
        modified to include metadata.

        :param str fn: The full path to the file.

        :return: None.
        """
        # Store data as dictionaries.
        bodies = {}
        indx = 0

        # Dictionaries to attach wing reference surfaces to wing bodies using
        # reference surface ID as the key.
        wing_bodies = {}
        ref_surfs = {}

        # Build a compound for geometric sets.
        compound = TopoDS_Compound()
        BRep_Builder().MakeCompound(compound)

        # Read file with OCCT.
        step_reader = STEPControl_Reader()
        status = step_reader.ReadFile(fn)
        if status not in [IFSelect_RetVoid, IFSelect_RetDone]:
            return bodies

        # Convert to desired units.
        Interface_Static.SetCVal_("xstep.cascade.unit", Settings.units)

        # Check.
        failsonly = False
        step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)

        # Transfer. OpenVSP STEP files result in one root and one shape (a
        # compound).
        step_reader.TransferRoot(1)
        master_shape = step_reader.Shape(1)

        # Things needed to extract names from STEP entities.
        session = step_reader.WS()
        transfer_reader = session.TransferReader()

        # Iterate over master shape to find compounds for geometric sets. These
        # sets contain the metadata and the surfaces that make up the
        # component.
        iterator = TopoDS_Iterator(master_shape, True, True)
        more = True
        while iterator.More() and more:
            # The compound.
            compound = iterator.Value()
            # Hack to handle single component for now...
            if compound.ShapeType() != TopAbs_COMPOUND:
                compound = master_shape
                more = False
            # Get the metadata.
            rep_item = transfer_reader.EntityFromShapeResult(compound, 1)
            name = rep_item.Name().ToCString()

            # Unnamed body
            if not name:
                indx += 1
                comp_name = '.'.join(['Body', str(indx)])
                msg = ' '.join(['---Processing OpenVSP component:', comp_name])
                logger.info(msg)
                solid = _build_solid(compound, self._divide)
                if solid:
                    body = Body(solid)
                    bodies[comp_name] = body
                iterator.Next()
                continue
            metadata = json.loads(name)

            # Process wing reference surface and continue.
            key = 'm_SurfType'
            if key in metadata and metadata[key] == 99:
                # Get surface
                sref = ImportVSP.process_sref(compound)
                # Get Sref ID
                sref_id = metadata['ID']
                ref_surfs[sref_id] = sref
                # Next shape.
                iterator.Next()
                continue

            comp_name = metadata['m_Name']
            if comp_name in bodies:
                indx += 1
                comp_name = '.'.join([comp_name, str(indx)])

            # Process component.
            msg = ' '.join(['---Processing OpenVSP component:', comp_name])
            logger.info(msg)

            # Wing
            if metadata['m_Type'] == 5 and metadata['m_SurfType'] != 99:
                wing = _process_wing(compound, self._divide, self._restrict,
                                     self._tol, self._reloft)
                if wing is not None:
                    wing.set_label(comp_name)
                    bodies[comp_name] = wing
                    sref_id = metadata['Sref ID']
                    wing_bodies[sref_id] = wing

            # Fuselage
            elif metadata['m_Type'] in [4, 9]:
                fuse = _process_fuse(compound, self._divide)
                if fuse is not None:
                    fuse.set_label(comp_name)
                    bodies[comp_name] = fuse

            # Unknown
            else:
                solid = _build_solid(compound, self._divide)
                if solid:
                    body = Body(solid)
                    bodies[comp_name] = body

            # Next shape.
            iterator.Next()

        # Attach wing reference surfaces to the bodies.
        for sref_id in wing_bodies:
            if sref_id not in ref_surfs:
                continue
            wing = wing_bodies[sref_id]
            sref = ref_surfs[sref_id]
            wing.set_sref(sref)

        # Update
        self._bodies.update(bodies)

    @staticmethod
    def rebuild_wing_solid(srfs, divide_closed=True, reloft=False, tol=0.01):
        """
        Rebuild a solid shape from the OpenVSP wing surface(s). If only one
        surface is provided then it is assumed that a single surface models the
        OML and it will be split and modified at the root, tip, and trailing
        edge. This single surface should have similar form and
        parametrization as the original OpenVSP surface. If more than once
        surface is provided then it is assumed that the surfaces were split
        during OpenVSP export and are simply sewn together to form the solid.

        :param srfs: The wing surface(s) used to rebuild the solid.
        :type srfs: collections.Sequence(afem.geometry.entities.Surface)
        :param bool divide_closed: Option to divide closed faces.
        :param bool reloft: For wings that are not split, this option will
            extract isocurves at each spanwise cross section, tessellate the
            curve using *tol*, and then approximate the section using a C1
            continuous curve. These curves are then used to generate a new
            wing surface. This method is experimental.
        :param float tol: Tolerance for approximation if
            *bspline_restrict* or *reloft* is *True*.

        :return: The new solid.
        :rtype: OCCT.TopoDS.TopoDS_Solid

        :raise ValueError: If no surfaces are provided.
        """
        faces = [FaceBySurface(s).face for s in srfs]
        compound = CompoundByShapes(faces).compound

        nsrfs = len(srfs)
        if nsrfs == 1:
            return _process_unsplit_wing(compound, divide_closed, reloft, tol)
        elif nsrfs > 1:
            return _build_solid(compound, divide_closed)
        else:
            raise ValueError('No surfaces provided.')

    @staticmethod
    def process_sref(compound):
        """
        Process a wing reference surface. The compound should contain a single
        face. The underlying surface of this faces will be used to generate a
        new NurbsSurface. Since the OpenVSP surface uses uniform
        parametrization, a chord line is extracted at each unique knot in the
        v-direction. Then a linear surface is lofted through these lines to
        form a bilinear surface where the u-direction is chordwise and the
        v-direction is spanwise.

        :param OCCT.TopoDS.TopoDS_Compound compound: The compound.

        :return: The reference surface.
        :rtype: afem.geometry.entities.NurbsSurface

        :raise TypeError: If the underlying surface cannot be downcasted to a
            NurbsSurface.
        """
        # Get underlying surface
        top_exp = TopExp_Explorer(compound, TopAbs_FACE)
        face = CheckShape.to_face(top_exp.Current())
        hsrf = BRep_Tool.Surface_(face)

        # Create NurbsSurface
        srf = NurbsSurface(hsrf)

        # Loft new surface
        vknots = srf.vknots
        crvs = []
        for vi in vknots:
            c = srf.v_iso(vi)
            crvs.append(c)
        srf = NurbsSurfaceByInterp(crvs, 1).surface
        return srf

    @staticmethod
    def rebuild_wing_sref(srf):
        """
        Attempt to rebuild a wing reference surface using the given OML
        surface. This method is intended to operate on OpenVSP surfaces that
        define the OML of a wing component and/or have similar parametrization.

        :param afem.geometry.entities.NurbsSurface srf: The surface.

        :return: The reference surface.
        :rtype: afem.geometry.entities.NurbsSurface
        """
        # For VSP, wing surfaces are degree=1 in the spanwise direction, so
        # each knot vector usually represents where a wing cross section is
        # located. For each knot vector in spanwise direction, generate a
        # straight line segment between the LE and TE (i.e., the chord line).
        # These chords will serve as the wing reference surface. This assumes
        # that the LE is at v=0.5 in the local parametric domain.
        uknots = srf.uknots
        le_param = srf.local_to_global_param('v', 0.5)
        te_param = srf.local_to_global_param('v', 0.)
        chords = []
        for u in uknots:
            le = srf.eval(u, le_param)
            te = srf.eval(u, te_param)
            c = NurbsCurveByPoints([le, te]).curve
            chords.append(c)

        # VSP wing components wrap around the ends and there may be duplicate
        # chord lines at the root and tip. Retain only unique lines based on LE
        # point.
        unique_chords = [chords[0]]
        for i in range(1, len(chords)):
            c0 = unique_chords[-1]
            ci = chords[i]
            if not ci.eval(0.).is_equal(c0.eval(0.)):
                unique_chords.append(ci)

        # Create degree=1 reference surface by skinning chord lines
        sref = NurbsSurfaceByInterp(unique_chords, 1).surface

        # Set domains to be 0 to 1.
        sref.set_udomain(0., 1.)
        sref.set_vdomain(0., 1.)

        return sref


def _build_solid(compound, divide_closed):
    # Get all the faces in the compound. The surfaces must be split. Discard
    # any with zero area.
    top_exp = TopExp_Explorer(compound, TopAbs_FACE)
    faces = []
    while top_exp.More():
        shape = top_exp.Current()
        face = CheckShape.to_face(shape)
        fprop = GProp_GProps()
        BRepGProp.SurfaceProperties_(face, fprop, 1.0e-7)
        a = fprop.Mass()
        if a <= 1.0e-7:
            top_exp.Next()
            continue
        faces.append(face)
        top_exp.Next()

    # Replace any planar B-Spline surfaces with planes.
    non_planar_faces = []
    planar_faces = []
    for f in faces:
        hsrf = BRep_Tool.Surface_(f)
        try:
            is_pln = GeomLib_IsPlanarSurface(hsrf, 1.0e-7)
            if is_pln.IsPlanar():
                w = ShapeAnalysis.OuterWire_(f)
                # Fix the wire because they are usually degenerate edges in
                # the planar end caps.
                builder = BRepBuilderAPI_MakeWire()
                for e in ExploreShape.get_edges(w):
                    if LinearProps(e).length > 1.0e-7:
                        builder.Add(e)
                w = builder.Wire()
                fix = ShapeFix_Wire()
                fix.Load(w)
                geom_pln = Geom_Plane(is_pln.Plan())
                fix.SetSurface(geom_pln)
                fix.FixReorder()
                fix.FixConnected()
                fix.FixEdgeCurves()
                fix.FixDegenerated()
                w = fix.WireAPIMake()
                # Build the planar face.
                fnew = BRepBuilderAPI_MakeFace(w, True).Face()
                planar_faces.append(fnew)
            else:
                non_planar_faces.append(f)
        except RuntimeError:
            logger.info('Failed to check for planar face...')
            non_planar_faces.append(f)

    # Make a compound of the faces
    shape = CompoundByShapes(non_planar_faces + planar_faces).compound

    # Split closed faces
    if divide_closed:
        divide = ShapeUpgrade_ShapeDivideClosed(shape)
        divide.Perform()
        shape = divide.Result()

    # Sew shape
    sew = BRepBuilderAPI_Sewing(1.0e-7)
    sew.Load(shape)
    sew.Perform()
    sewn_shape = sew.SewedShape()

    if sewn_shape.ShapeType() == TopAbs_FACE:
        face = sewn_shape
        sewn_shape = TopoDS_Shell()
        builder = BRep_Builder()
        builder.MakeShell(sewn_shape)
        builder.Add(sewn_shape, face)

    # Attempt to unify planar domains
    unify_shp = ShapeUpgrade_UnifySameDomain(sewn_shape, False, False, False)
    try:
        unify_shp.UnifyFaces()
        shape = unify_shp.Shape()
    except RuntimeError:
        logger.info('...failed to unify faces on solid.')
        shape = sewn_shape

    # Make solid
    shell = ExploreShape.get_shells(shape)[0]
    solid = ShapeFix_Solid().SolidFromShell(shell)

    # Limit tolerance
    FixShape.limit_tolerance(solid)

    # Check the solid and attempt to fix
    check_shp = BRepCheck_Analyzer(solid, True)
    if not check_shp.IsValid():
        logger.info('\tFixing the solid...')
        fix = ShapeFix_Solid(solid)
        fix.Perform()
        solid = fix.Solid()
        check_shp = BRepCheck_Analyzer(solid, True)
        if not check_shp.IsValid():
            logger.info('\t...solid could not be fixed.')
            logger.info('\tShape diagnostics:')
            _topods_iterator_check(solid, check_shp)
    else:
        tol = ExploreShape.global_tolerance(solid)
        logger.info(
            '\tSuccessfully generated solid with tolerance={}'.format(tol))

    return solid


def _process_wing(compound, divide_closed, bspline_restrict, tol, reloft):
    # Note that for VSP wings, the spanwise direction is u and the chord
    # direction is v, where v=0 is the TE and follows the lower surface fwd to
    # the LE, and then aft along the upper surface to the TE.

    # Process based on number of faces in compound assuming split/no split
    # option was used.
    faces = ExploreShape.get_faces(compound)
    vsp_surf = None
    if len(faces) == 1:
        solid = _process_unsplit_wing(compound, divide_closed, reloft, tol)
        vsp_surf = ExploreShape.surface_of_face(faces[0])
    else:
        solid = _build_solid(compound, divide_closed)

    if not solid:
        return None

    if bspline_restrict:
        solid = _bspline_restrict(solid, tol)

    wing = Body(solid)

    if vsp_surf:
        vsp_surf = NurbsSurface(vsp_surf.object)
        wing.add_metadata('vsp surface', vsp_surf)
        upr_srf = vsp_surf.copy()
        v_le = vsp_surf.local_to_global_param('v', 0.5)
        upr_srf.segment(vsp_surf.u1, vsp_surf.u2, v_le, vsp_surf.v2)
        wing.add_metadata('upper surface', upr_srf)
        lwr_srf = vsp_surf.copy()
        lwr_srf.segment(vsp_surf.u1, vsp_surf.u2, vsp_surf.v1, v_le)
        wing.add_metadata('lower surface', lwr_srf)

    return wing


def _process_fuse(compound, divide_closed):
    # For VSP fuselages, the longitudinal direction is u, and the
    # circumferential direction is v.
    # Build the solid.
    solid = _build_solid(compound, divide_closed)
    if not solid:
        return None

    fuselage = Body(solid)

    faces = ExploreShape.get_faces(compound)
    if len(faces) == 1:
        vsp_surf = ExploreShape.surface_of_face(faces[0])
        vsp_surf = NurbsSurface(vsp_surf.object)
        fuselage.add_metadata('vsp surface', vsp_surf)

    return fuselage


def _process_unsplit_wing(compound, divide_closed, reloft, tol):
    # Process a wing that was generated without "Split Surfs" option.

    faces = ExploreShape.get_faces(compound)
    if len(faces) != 1:
        return None
    face = faces[0]

    # Get the surface.
    master_surf = ExploreShape.surface_of_face(face)
    master_surf = NurbsSurface(master_surf.object)
    uknots, vknots = master_surf.uknots, master_surf.vknots
    vsplit = master_surf.local_to_global_param('v', 0.5)

    # Segment off the end caps and the trailing edges.
    u1, u2 = uknots[1], uknots[-2]
    v1, v2 = vknots[1], vknots[-2]
    s1 = master_surf.copy()
    s1.segment(u1, u2, v1, v2)

    # Reloft the surface by tessellating a curve at each spanwise knot. This
    # enforces C1 continuity but assumes linear spanwise wing which may not
    # support blending wing sections in newer versions of OpenVSP. Also, since
    # the tessellated curves may not match up to the wing end caps making sewing
    # unreliable, flat end caps are assumed.
    if reloft:
        s1 = _reloft_wing_surface(s1, tol)

        # Generate new flat end caps using isocurves at the root and tip of this
        # new surface
        c0 = s1.v_iso(s1.v1)
        c1 = s1.v_iso(s1.v2)
        e0 = EdgeByCurve(c0).edge
        e1 = EdgeByCurve(c1).edge
        w0 = WireByEdges(e0).wire
        w1 = WireByEdges(e1).wire

        f0 = BRepBuilderAPI_MakeFace(w0, True).Face()
        f1 = BRepBuilderAPI_MakeFace(w1, True).Face()

        # Make faces of surfaces
        f = FaceBySurface(s1).face
        new_faces = [f, f0, f1]
    else:
        # Reparamterize knots in spanwise direction to be chord length instead
        # of uniform. Use isocurve at quarter-chord to determine knot values.
        # This only works as long as surfaces are linear.
        c0 = NurbsCurve.downcast(s1.u_iso(s1.u1))
        c0.segment(vsplit, c0.u2)
        qc_u = PointFromParameter(c0, vsplit, 0.25 * c0.length).parameter
        c = NurbsCurve.downcast(s1.v_iso(qc_u))
        pnts = [c.eval(u) for u in c.knots]
        new_uknots = chord_parameters(pnts, 0., 1.)
        s1.set_uknots(new_uknots)

        # Segment off end caps
        u1, u2 = uknots[0], uknots[1]
        v1, v2 = vknots[1], vsplit
        s2 = master_surf.copy()
        s2.segment(u1, u2, v1, v2)

        u1, u2 = uknots[0], uknots[1]
        v1, v2 = vsplit, vknots[-2]
        s3 = master_surf.copy()
        s3.segment(u1, u2, v1, v2)

        u1, u2 = uknots[-2], uknots[-1]
        v1, v2 = vknots[1], vsplit
        s4 = master_surf.copy()
        s4.segment(u1, u2, v1, v2)

        u1, u2 = uknots[-2], uknots[-1]
        v1, v2 = vsplit, vknots[-2]
        s5 = master_surf.copy()
        s5.segment(u1, u2, v1, v2)

        # Make faces of surfaces
        new_faces = []
        for s in [s1, s2, s3, s4, s5]:
            f = BRepBuilderAPI_MakeFace(s.object, 0.).Face()
            new_faces.append(f)

    # Segment off TE.
    u1, u2 = uknots[0], uknots[-1]
    v1, v2 = vknots[0], vknots[1]
    s6 = master_surf.copy()
    s6.segment(u1, u2, v1, v2)

    u1, u2 = uknots[0], uknots[-1]
    v1, v2 = vknots[-2], vknots[-1]
    s7 = master_surf.copy()
    s7.segment(u1, u2, v1, v2)

    # Split the TE surface at each u-knot.
    usplits = TColStd_HSequenceOfReal()
    for ui in uknots[1:-1]:
        usplits.Append(ui)

    split = ShapeUpgrade_SplitSurface()
    split.Init(s6.object)
    split.SetUSplitValues(usplits)
    split.Perform()
    comp_surf1 = split.ResSurfaces()

    split = ShapeUpgrade_SplitSurface()
    split.Init(s7.object)
    split.SetUSplitValues(usplits)
    split.Perform()
    comp_surf2 = split.ResSurfaces()

    # For each patch in the composite surfaces create a face.
    for i in range(1, comp_surf1.NbUPatches() + 1):
        for j in range(1, comp_surf1.NbVPatches() + 1):
            hpatch = comp_surf1.Patch(i, j)
            f = BRepBuilderAPI_MakeFace(hpatch, 0.).Face()
            new_faces.append(f)

    for i in range(1, comp_surf2.NbUPatches() + 1):
        for j in range(1, comp_surf2.NbVPatches() + 1):
            hpatch = comp_surf2.Patch(i, j)
            f = BRepBuilderAPI_MakeFace(hpatch, 0.).Face()
            new_faces.append(f)

    # Put all faces into a compound a generate solid.
    new_compound = CompoundByShapes(new_faces).compound

    return _build_solid(new_compound, divide_closed)


def _bspline_restrict(solid, tol):
    """
    Attempt to re-fit OpenVSP surfaces with more continuity.
    """
    # Use dmax=1 because that was only way to get the tool actually refit the
    # surfaces other than another surface where the multiplicity equaled the
    # degree. Not sure why the tool operates this way.
    logger.info('\tApplying ShapeBSplineRestriction tool...')
    tool = ShapeBSplineRestriction(solid, dmax=1, tol3d=tol)
    if not tool.is_done:
        logger.info('Method unsuccessful. Using original solid.')
        return solid

    # Get new shape and solid
    new_shape = tool.modified_shape(solid)
    new_solid = CheckShape.to_solid(new_shape)

    # Limit/fix tolerance
    FixShape.limit_tolerance(new_solid)
    tol = ExploreShape.global_tolerance(new_solid)
    if not CheckShape.is_valid(new_solid):
        logger.info('Shape invalid. Using original solid.')
        return solid

    logger.info('\tMethod successful with surface error: {}'.format(
        tool.error_surface))
    logger.info('\tNew shape tolerance: {}'.format(tol))

    return new_solid


def _reloft_wing_surface(srf, tol):
    """
    Attempt to reloft an OpenVSP wing surface which was not split to achieve
    higher continuity.
    """
    logger.info('\tAttempting to reloft the surface...')
    # Gather isocurves at each section, tessellate, and approximate
    crvs = []
    for u in srf.uknots:
        c0 = srf.u_iso(u)
        adp_crv = GeomAdaptor_Curve(c0.object)
        tool = GCPnts_QuasiUniformDeflection(adp_crv, tol)
        if not tool.IsDone():
            logger.info('\tTessellation failed. Using original surface.')
            return srf
        pnts = [c0.eval(tool.Parameter(i)) for i in
                range(1, tool.NbPoints() + 1)]
        c = NurbsCurveByApprox(pnts, tol=tol, continuity='C1').curve
        crvs.append(c)
    return NurbsSurfaceByInterp(crvs, 1).surface


def _topods_iterator_check(shape, check):
    """
    Iterate on the shape and dump errors.
    """
    it = TopoDS_Iterator(shape)
    while it.More():
        sub_shape = it.Value()
        result = check.Result(sub_shape)
        list_of_status = result.Status()
        for status in list_of_status:
            if status != BRepCheck_NoError:
                str_log = '\t\t{0}-->{1}'.format(status, sub_shape.ShapeType())
                logger.info(str_log)
        it.Next()
        _topods_iterator_check(sub_shape, check)
