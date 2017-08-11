Topology
========
This section describes the topology package. The topology tools can be
imported by::

    from afem.topology import *

Entities
--------
.. py:currentmodule:: OCC.TopoDS

OpenCASCADE entities are provided here for reference. Rarely will the user
ever need to initialize any of these entities. The topology builders should
be used instead.

Bounding Box
~~~~~~~~~~~~
.. autoclass:: afem.topology.entities.BBox

TopoDS_Shape
~~~~~~~~~~~~
.. autoclass:: TopoDS_Shape

TopoDS_Vertex
~~~~~~~~~~~~~
.. autoclass:: TopoDS_Vertex

TopoDS_Edge
~~~~~~~~~~~
.. autoclass:: TopoDS_Edge

TopoDS_Wire
~~~~~~~~~~~
.. autoclass:: TopoDS_Wire

TopoDS_Face
~~~~~~~~~~~
.. autoclass:: TopoDS_Face

TopoDS_Shell
~~~~~~~~~~~~
.. autoclass:: TopoDS_Shell

TopoDS_Solid
~~~~~~~~~~~~
.. autoclass:: TopoDS_Solid

TopoDS_CompSolid
~~~~~~~~~~~~~~~~
.. autoclass:: TopoDS_CompSolid

TopoDS_Compound
~~~~~~~~~~~~~~~
.. autoclass:: TopoDS_Compound

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

WireBySplit
~~~~~~~~~~~
.. autoclass:: WireBySplit

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

ExploreShape
~~~~~~~~~~~~
.. autoclass:: ExploreShape

ExploreWire
~~~~~~~~~~~
.. autoclass:: ExploreWire

ExploreFreeEdges
~~~~~~~~~~~~~~~~
.. autoclass:: ExploreFreeEdges

Modify
------
.. py:currentmodule:: afem.topology.modify

FixShape
~~~~~~~~
.. autoclass:: FixShape

DivideClosedShape
~~~~~~~~~~~~~~~~~
.. autoclass:: DivideClosedShape

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

Boolean
-------
.. py:currentmodule:: afem.topology.bop

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
