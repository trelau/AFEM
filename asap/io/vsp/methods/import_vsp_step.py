from __future__ import print_function

import json

import OCC.ShapeAnalysis as ShapeAnalysis
from OCC.BRep import BRep_Builder, BRep_Tool
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_Sewing
from OCC.BRepCheck import BRepCheck_Analyzer
from OCC.BRepGProp import brepgprop
from OCC.GProp import GProp_GProps
from OCC.GeomAdaptor import GeomAdaptor_Surface
from OCC.GeomLib import GeomLib_IsPlanarSurface
from OCC.IFSelect import IFSelect_ItemsByEntity
from OCC.Interface import Interface_Static
from OCC.STEPControl import STEPControl_Reader
from OCC.ShapeFix import ShapeFix_Shape, ShapeFix_Solid
from OCC.ShapeUpgrade import ShapeUpgrade_ShapeDivideClosed, \
    ShapeUpgrade_SplitSurface, ShapeUpgrade_UnifySameDomain
from OCC.StepRepr import StepRepr_RepresentationItem
from OCC.TColStd import TColStd_HSequenceOfReal
from OCC.TopAbs import TopAbs_FACE
from OCC.TopExp import TopExp_Explorer
from OCC.TopoDS import TopoDS_Compound, TopoDS_Iterator, TopoDS_Shell

from ....geometry import CreateGeom
from ....geometry.methods.create import create_nurbs_surface_from_occ
from ....oml.body import Body
from ....oml.fuselage import Fuselage
from ....oml.wing import Wing
from ....topology import ShapeTools


def import_vsp_step(fname, divide_closed):
    """
    Import a OpenVSP STEP file and translate into ASAP wing and fuselage
    bodies.

    :param str fname: Filename (and path) to STEP file.
    :param divide_closed:

    :return: Dictionaries of wing and fuselage bodies.
    :rtype: dict
    """
    # Store data as dictionaries.
    bodies = {}
    indx = 0

    # Dictionaries to attach wing reference surfaces to wing bodies using
    # reference surface ID as the key.
    wing_bodies = {}
    ref_surfs = {}

    # Convert to millimeters.
    Interface_Static.SetCVal("xstep.cascade.unit", "MM")
    Interface_Static.SetIVal("read.step.ideas", 1)
    Interface_Static.SetIVal("read.step.nonmanifold", 1)

    # Build a compound for geometric sets.
    compound = TopoDS_Compound()
    BRep_Builder().MakeCompound(compound)

    # Read file with OCC.
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(fname)
    if status > 1:
        return bodies

    # Convert to inches.
    Interface_Static.SetCVal("xstep.cascade.unit", "INCH")

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
    transfer_reader = session.GetObject().TransferReader().GetObject()

    # Iterate over master shape to find compounds for geometric sets. These
    # sets contain the metadata and the surfaces that make up the component.
    # TODO How to handle a single component?
    iterator = TopoDS_Iterator(master_shape, True, True)
    while iterator.More():
        # The compound.
        compound = iterator.Value()
        # Get the metadata.
        entity = transfer_reader.EntityFromShapeResult(compound, 1)
        rep_item = StepRepr_RepresentationItem().GetHandle().DownCast(entity)
        name = rep_item.GetObject().Name().GetObject().ToCString()
        if not name:
            indx += 1
            comp_name = '.'.join(['Body', str(indx)])
            print('---Processing OpenVSP component:', comp_name)
            solid = _build_solid(compound, divide_closed)
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
            sref = _process_sref(compound)
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
        print('---Processing OpenVSP component:', comp_name)

        # Wing
        if metadata['m_Type'] == 5 and metadata['m_SurfType'] != 99:
            wing = _process_wing(compound, divide_closed)
            if wing is not None:
                wing.set_name(comp_name)
                bodies[comp_name] = wing
                sref_id = metadata['Sref ID']
                wing_bodies[sref_id] = wing

        # Fuselage
        elif metadata['m_Type'] in [4, 9]:
            fuse = _process_fuse(compound, divide_closed)
            if fuse is not None:
                fuse.set_name(comp_name)
                bodies[comp_name] = fuse

        # Next shape.
        iterator.Next()

    # Attach wing reference surfaces to the bodies.
    for sref_id in wing_bodies:
        if sref_id not in ref_surfs:
            continue
        wing = wing_bodies[sref_id]
        sref = ref_surfs[sref_id]
        status = wing.set_sref(sref)
        if not status:
            print('Failed to set wing reference surface for a wing body.')

    return bodies


def _build_solid(compound, divide_closed):
    # Get all the faces in the compound. The surfaces must be split. Discard
    # any with zero area.
    top_exp = TopExp_Explorer(compound, TopAbs_FACE)
    faces = []
    while top_exp.More():
        shape = top_exp.Current()
        face = ShapeTools.to_face(shape)
        fprop = GProp_GProps()
        brepgprop.SurfaceProperties(face, fprop, 1.0e-7)
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
        hsrf = BRep_Tool.Surface(f)
        try:
            is_pln = GeomLib_IsPlanarSurface(hsrf, 1.0e-7)
            if is_pln.IsPlanar():
                w = ShapeAnalysis.shapeanalysis_OuterWire(f)
                # Fix the wire because they are usually degenerate edges in
                # the planar end caps.
                fix = ShapeFix_Shape(w)
                fix.Perform()
                w = ShapeTools.to_wire(fix.Shape())
                # Build the planar face.
                fnew = BRepBuilderAPI_MakeFace(is_pln.Plan(), w, True).Face()
                planar_faces.append(fnew)
            else:
                non_planar_faces.append(f)
        except RuntimeError:
            print('Failed to check for planar face...')
            non_planar_faces.append(f)

    # # Check the faces.
    # for i, f in enumerate(non_planar_faces):
    #     check = BRepCheck_Analyzer(f, True)
    #     if not check.IsValid():
    #         print('Non-planar face is not valid...')
    #         fix = ShapeFix_Face(f)
    #         fix.Perform()
    #         fnew = fix.Result()
    #         check = BRepCheck_Analyzer(fnew, True)
    #         if not check.IsValid():
    #             print('...face could not be fixed.')
    #         else:
    #             non_planar_faces[i] = fnew
    #
    # for i, f in enumerate(planar_faces):
    #     check = BRepCheck_Analyzer(f, True)
    #     if not check.IsValid():
    #         print('Planar face is not valid...')
    #         fix = ShapeFix_Face(f)
    #         fix.Perform()
    #         fnew = fix.Result()
    #         check = BRepCheck_Analyzer(fnew, True)
    #         if not check.IsValid():
    #             print('...face could not be fixed.')
    #         else:
    #             planar_faces[i] = fnew

    # Make a shell and a solid.
    shell = TopoDS_Shell()
    builder = BRep_Builder()
    builder.MakeShell(shell)
    for f in non_planar_faces + planar_faces:
        builder.Add(shell, f)

    # Sew shape.
    sew = BRepBuilderAPI_Sewing(0.01)
    sew.Load(shell)
    sew.Perform()
    sewn_shape = sew.SewedShape()

    if sewn_shape.ShapeType() == TopAbs_FACE:
        face = sewn_shape
        sewn_shape = TopoDS_Shell()
        builder = BRep_Builder()
        builder.MakeShell(sewn_shape)
        builder.Add(sewn_shape, face)

    # Attempt to unify planar domains.
    unify_shp = ShapeUpgrade_UnifySameDomain(sewn_shape, False, False, False)
    try:
        unify_shp.UnifyFaces()
        shape = unify_shp.Shape()
    except RuntimeError:
        shape = sewn_shape

    # Make solid.
    shell = ShapeTools.to_shell(shape)
    solid = ShapeFix_Solid().SolidFromShell(shell)

    # Split closed faces of the solid to make OCC more robust.
    if divide_closed:
        divide = ShapeUpgrade_ShapeDivideClosed(solid)
        divide.Perform()
        solid = divide.Result()

    # Make sure it's a solid.
    solid = ShapeTools.to_solid(solid)

    # Check the solid and attempt to fix.
    check_shp = BRepCheck_Analyzer(solid, True)
    if not check_shp.IsValid():
        print('Fixing the solid...')
        fix = ShapeFix_Solid(solid)
        fix.Perform()
        solid = fix.Solid()
        check_shp = BRepCheck_Analyzer(solid, True)
        if not check_shp.IsValid():
            print('...solid could not be fixed.')
    else:
        print('Successfully generated solid.')

    return solid


def _process_wing(compound, divide_closed):
    # Note that for VSP wings, the spanwise direction is u and the chord
    # direction is v, where v=0 is the TE and follows the lower surface fwd to
    # the LE, and then aft along the upper surface to the TE.

    # Process based on number of faces in compound assuming split/no split
    # option was used.
    faces = ShapeTools.get_faces(compound)
    vsp_surf = None
    if len(faces) == 1:
        solid = _process_unsplit_wing(compound, divide_closed)
        vsp_surf = ShapeTools.surface_of_face(faces[0])
    else:
        solid = _build_solid(compound, divide_closed)

    if not solid:
        return None

    wing = Wing(solid)

    if vsp_surf:
        wing.add_metadata('vsp surface', vsp_surf)
        upr_srf = CreateGeom.copy_geom(vsp_surf)
        v_le = vsp_surf.local_to_global_param('v', 0.5)
        upr_srf.segment(vsp_surf.u1, vsp_surf.u2, v_le, vsp_surf.v2)
        wing.add_metadata('upper surface', upr_srf)
        lwr_srf = CreateGeom.copy_geom(vsp_surf)
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

    fuselage = Fuselage(solid)

    faces = ShapeTools.get_faces(compound)
    if len(faces) == 1:
        vsp_surf = ShapeTools.surface_of_face(faces[0])
        fuselage.add_metadata('vsp surface', vsp_surf)

    return fuselage


def _process_sref(compound):
    """
    Process a wing reference surface.
    """
    # Wing reference surfaces should be a single bilinear surface. Extract
    # this from the compound.
    top_exp = TopExp_Explorer(compound, TopAbs_FACE)
    face = ShapeTools.to_face(top_exp.Current())
    hsrf = BRep_Tool.Surface(face)
    adp_srf = GeomAdaptor_Surface(hsrf)
    occ_srf = adp_srf.BSpline().GetObject()

    # Convert to ASAP NurbsSurface.
    srf = create_nurbs_surface_from_occ(occ_srf)

    # The surface originally has uniform parameterization from OpenVSP,
    # extract curves at each knot value and wing_skin a C0 surface using chord
    # length parameterization.
    vknots = srf.vknots
    crvs = []
    for vi in vknots:
        c = CreateGeom.isocurve(srf, None, vi)
        crvs.append(c)
    srf = CreateGeom.interp_curves_by_surface(crvs, 1, 'chord')
    return srf


def _process_unsplit_wing(compound, divide_closed):
    # Process a wing that was generated without "Split Surfs" option.

    faces = ShapeTools.get_faces(compound)
    if len(faces) != 1:
        return None
    face = faces[0]

    # Get the surface.
    master_surf = ShapeTools.surface_of_face(face)
    uknots, vknots = master_surf.uknots, master_surf.vknots
    vsplit = master_surf.local_to_global_param('v', 0.5)

    # Segment off the end caps and the trailing edges.
    u1, u2 = uknots[1], uknots[-2]
    v1, v2 = vknots[1], vknots[-2]
    s1 = CreateGeom.copy_geom(master_surf)
    s1.segment(u1, u2, v1, v2)

    # Segment off end caps and the trailing edge and split at LE.
    u1, u2 = uknots[0], uknots[1]
    v1, v2 = vknots[1], vsplit
    s2 = CreateGeom.copy_geom(master_surf)
    s2.segment(u1, u2, v1, v2)

    u1, u2 = uknots[0], uknots[1]
    v1, v2 = vsplit, vknots[-2]
    s3 = CreateGeom.copy_geom(master_surf)
    s3.segment(u1, u2, v1, v2)

    u1, u2 = uknots[-2], uknots[-1]
    v1, v2 = vknots[1], vsplit
    s4 = CreateGeom.copy_geom(master_surf)
    s4.segment(u1, u2, v1, v2)

    u1, u2 = uknots[-2], uknots[-1]
    v1, v2 = vsplit, vknots[-2]
    s5 = CreateGeom.copy_geom(master_surf)
    s5.segment(u1, u2, v1, v2)

    # Make faces of surface.
    new_faces = []
    for s in [s1, s2, s3, s4, s5]:
        f = BRepBuilderAPI_MakeFace(s.GetHandle(), 0.).Face()
        new_faces.append(f)

    # Segment off TE.
    u1, u2 = uknots[0], uknots[-1]
    v1, v2 = vknots[0], vknots[1]
    s6 = CreateGeom.copy_geom(master_surf)
    s6.segment(u1, u2, v1, v2)

    u1, u2 = uknots[0], uknots[-1]
    v1, v2 = vknots[-2], vknots[-1]
    s7 = CreateGeom.copy_geom(master_surf)
    s7.segment(u1, u2, v1, v2)

    # Split the TE surface at each uknot.
    usplits = TColStd_HSequenceOfReal()
    for ui in uknots[1:-1]:
        usplits.Append(ui)

    split = ShapeUpgrade_SplitSurface()
    split.Init(s6.GetHandle())
    split.SetUSplitValues(usplits.GetHandle())
    split.Perform()
    comp_surf1 = split.ResSurfaces().GetObject()

    split = ShapeUpgrade_SplitSurface()
    split.Init(s7.GetHandle())
    split.SetUSplitValues(usplits.GetHandle())
    split.Perform()
    comp_surf2 = split.ResSurfaces().GetObject()

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
    new_compound = ShapeTools.make_compound(new_faces)

    return _build_solid(new_compound, divide_closed)
