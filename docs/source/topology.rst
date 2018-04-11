Topology
========
The ``topology`` package provides tools for the creation and use of OpenCASCADE
native topology (i.e., shapes). While geometry defines curves and surfaces,
topology describes their connectivity and boundary representation. OpenCASCADE
shapes are the core building blocks for building more complex parts and
assemblies. The topology entities and tools can be imported by::

    from afem.topology import *

It is important to note that AFEM does not provide a "Pythonic" wrapper layer
for OpenCASCADE shapes like it does for geometry. The user should become
familiar with OpenCASCADE topology by reviewing OpenCASCADE documentation:

* `OpenCASCADE Reference Manual-Topology <https://www.opencascade.com/doc/occt-7.2.0/overview/html/occt_user_guides__modeling_data.html#occt_modat_5>`_

Many tools are provided for user to easily create, modify, and operate on
OpenCASCADE shapes. If a type or tool is not available via AFEM, the user can
still use the **pyOCCT** package required by AFEM.

Entities
--------
.. py:currentmodule:: afem.topology.entities

Shape
~~~~~
.. autoclass:: Shape

Vertex
~~~~~~~
.. autoclass:: Vertex

Edge
~~~~
.. autoclass:: Edge

Wire
~~~~
.. autoclass:: Wire

Face
~~~~
.. autoclass:: Face

Shell
~~~~~
.. autoclass:: Shell

Solid
~~~~~
.. autoclass:: Solid

CompSolid
~~~~~~~~~
.. autoclass:: CompSolid

Compound
~~~~~~~~
.. autoclass:: Compound

Bounding Box
~~~~~~~~~~~~
.. autoclass:: BBox

Create
------
.. py:currentmodule:: afem.topology.create

VertexByPoint
~~~~~~~~~~~~~
.. autoclass:: VertexByPoint

EdgeByPoints
~~~~~~~~~~~~
.. autoclass:: EdgeByPoints

EdgeByVertices
~~~~~~~~~~~~~~
.. autoclass:: EdgeByVertices

EdgeByCurve
~~~~~~~~~~~
.. autoclass:: EdgeByCurve

EdgeByDrag
~~~~~~~~~~
.. autoclass:: EdgeByDrag

EdgeByWireConcat
~~~~~~~~~~~~~~~~
.. autoclass:: EdgeByWireConcat

WireByEdges
~~~~~~~~~~~
.. autoclass:: WireByEdges

WiresByConnectedEdges
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: WiresByConnectedEdges

WireByPlanarOffset
~~~~~~~~~~~~~~~~~~
.. autoclass:: WireByPlanarOffset

WiresByShape
~~~~~~~~~~~~
.. autoclass:: WiresByShape

WireByPoints
~~~~~~~~~~~~~
.. autoclass:: WireByPoints

WireByConcat
~~~~~~~~~~~~
.. autoclass:: WireByConcat

FaceBySurface
~~~~~~~~~~~~~
.. autoclass:: FaceBySurface

FaceByPlane
~~~~~~~~~~~
.. autoclass:: FaceByPlane

FaceByPlanarWire
~~~~~~~~~~~~~~~~
.. autoclass:: FaceByPlanarWire

FaceByDrag
~~~~~~~~~~
.. autoclass:: FaceByDrag

ShellBySurface
~~~~~~~~~~~~~~
.. autoclass:: ShellBySurface

ShellByFaces
~~~~~~~~~~~~
.. autoclass:: ShellByFaces

ShellByDrag
~~~~~~~~~~~
.. autoclass:: ShellByDrag

ShellBySewing
~~~~~~~~~~~~~
.. autoclass:: ShellBySewing

SolidByShell
~~~~~~~~~~~~
.. autoclass:: SolidByShell

SolidByPlane
~~~~~~~~~~~~
.. autoclass:: SolidByPlane

SolidByDrag
~~~~~~~~~~~
.. autoclass:: SolidByDrag

CompoundByShapes
~~~~~~~~~~~~~~~~
.. autoclass:: CompoundByShapes

HalfspaceByShape
~~~~~~~~~~~~~~~~
.. autoclass:: HalfspaceByShape

HalfspaceBySurface
~~~~~~~~~~~~~~~~~~
.. autoclass:: HalfspaceBySurface

ShapeByFaces
~~~~~~~~~~~~
.. autoclass:: ShapeByFaces

ShapeByDrag
~~~~~~~~~~~
.. autoclass:: ShapeByDrag

BoxBuilder
~~~~~~~~~~
.. autoclass:: BoxBuilder

BoxBySize
~~~~~~~~~
.. autoclass:: BoxBySize

BoxBy2Points
~~~~~~~~~~~~
.. autoclass:: BoxBy2Points

CylinderByAxis
~~~~~~~~~~~~~~
.. autoclass:: CylinderByAxis

SphereByRadius
~~~~~~~~~~~~~~
.. autoclass:: SphereByRadius

SphereBy3Points
~~~~~~~~~~~~~~~
.. autoclass:: SphereBy3Points

PointAlongShape
~~~~~~~~~~~~~~~
.. autoclass:: PointAlongShape

PointsAlongShapeByNumber
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongShapeByNumber

PointsAlongShapeByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongShapeByDistance

PlaneByEdges
~~~~~~~~~~~~
.. autoclass:: PlaneByEdges

PlaneByIntersectingShapes
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlaneByIntersectingShapes

Explore
-------
.. py:currentmodule:: afem.topology.explore

ExploreWire
~~~~~~~~~~~
.. autoclass:: ExploreWire

ExploreFreeEdges
~~~~~~~~~~~~~~~~
.. autoclass:: ExploreFreeEdges

Modify
------
.. py:currentmodule:: afem.topology.modify

DivideClosedShape
~~~~~~~~~~~~~~~~~
.. autoclass:: DivideClosedShape

DivideContinuityShape
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DivideContinuityShape

DivideC0Shape
~~~~~~~~~~~~~~~~~
.. autoclass:: DivideC0Shape

UnifyShape
~~~~~~~~~~
.. autoclass:: UnifyShape

SewShape
~~~~~~~~
.. autoclass:: SewShape

RebuildShapeWithShapes
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RebuildShapeWithShapes

RebuildShapeByTool
~~~~~~~~~~~~~~~~~~
.. autoclass:: RebuildShapeByTool

RebuildShapesByTool
~~~~~~~~~~~~~~~~~~~
.. autoclass:: RebuildShapesByTool

ShapeBSplineRestriction
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: ShapeBSplineRestriction

Boolean
-------
.. py:currentmodule:: afem.topology.bop

BopCore
~~~~~~~
.. autoclass:: BopCore

BopAlgo
~~~~~~~
.. autoclass:: BopAlgo

FuseShapes
~~~~~~~~~~
.. autoclass:: FuseShapes

CutShapes
~~~~~~~~~
.. autoclass:: CutShapes

CommonShapes
~~~~~~~~~~~~
.. autoclass:: CommonShapes

IntersectShapes
~~~~~~~~~~~~~~~
.. autoclass:: IntersectShapes

SplitShapes
~~~~~~~~~~~~
.. autoclass:: SplitShapes

VolumesFromShapes
~~~~~~~~~~~~~~~~~
.. autoclass:: VolumesFromShapes

CutCylindricalHole
~~~~~~~~~~~~~~~~~~
.. autoclass:: CutCylindricalHole

LocalSplit
~~~~~~~~~~
.. autoclass:: LocalSplit

SplitShapeByEdges
~~~~~~~~~~~~~~~~~
.. autoclass:: SplitShapeByEdges

SplitWire
~~~~~~~~~~~
.. autoclass:: SplitWire

TrimOpenWire
~~~~~~~~~~~~
.. autoclass:: TrimOpenWire

Offset
------
.. py:currentmodule:: afem.topology.offset

ProjectShape
~~~~~~~~~~~~
.. autoclass:: ProjectShape

OffsetShape
~~~~~~~~~~~
.. autoclass:: OffsetShape

LoftShape
~~~~~~~~~
.. autoclass:: LoftShape

SweepShape
~~~~~~~~~~
.. autoclass:: SweepShape

SweepShapeWithNormal
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SweepShapeWithNormal

Distance
--------
.. py:currentmodule:: afem.topology.distance

DistanceShapeToShape
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistanceShapeToShape

DistanceShapeToShapes
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistanceShapeToShapes

Fix
---
.. py:currentmodule:: afem.topology.fix

FixShape
~~~~~~~~
.. autoclass:: FixShape

Properties
----------
.. py:currentmodule:: afem.topology.props

ShapeProps
~~~~~~~~~~
.. autoclass:: ShapeProps

LinearProps
~~~~~~~~~~~
.. autoclass:: LinearProps

SurfaceProps
~~~~~~~~~~~~
.. autoclass:: SurfaceProps

VolumeProps
~~~~~~~~~~~
.. autoclass:: VolumeProps

LengthOfShapes
~~~~~~~~~~~~~~
.. autoclass:: LengthOfShapes

AreaOfShapes
~~~~~~~~~~~~
.. autoclass:: AreaOfShapes

Check
-----
.. py:currentmodule:: afem.topology.check

CheckShape
~~~~~~~~~~
.. autoclass:: CheckShape

ClassifyPointInSolid
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: ClassifyPointInSolid

Transform
---------
.. automodule:: afem.topology.transform