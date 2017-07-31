Topology
========
This section describes the topology package. The topology tools can be
imported by::

    from afem.occ import *
    from afem.topology import *

Entities
--------
.. py:currentmodule:: OCC.TopoDS

OpenCASCADE entities are provided here for reference. Rarely will the user
ever need to initialize any of these entities. The topology builders should
be used instead.

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
.. py:currentmodule:: afem.occ.create

ShapeTools
~~~~~~~~~~
This class was never meant to be this big. It will be broken out in the near
future.

.. autoclass:: afem.topology.tools.ShapeTools

CreateShape
~~~~~~~~~~~
.. autoclass:: CreateShape

ShapeBuilder
~~~~~~~~~~~~
.. autoclass:: ShapeBuilder

CompoundByShapes
~~~~~~~~~~~~~~~~
.. autoclass:: CompoundByShapes

BoxByPlane
~~~~~~~~~~
.. autoclass:: BoxByPlane

WireByConnectedEdges
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: WireByConnectedEdges

ShapeBySweep
~~~~~~~~~~~~
.. autoclass:: ShapeBySweep

EdgeBySweep
~~~~~~~~~~~
.. autoclass:: EdgeBySweep

FaceBySweep
~~~~~~~~~~~
.. autoclass:: FaceBySweep

ShellBySweep
~~~~~~~~~~~~
.. autoclass:: ShellBySweep

SolidBySweep
~~~~~~~~~~~~
.. autoclass:: SolidBySweep

CompSolidBySweep
~~~~~~~~~~~~~~~~
.. autoclass:: CompSolidBySweep

WireByConcat
~~~~~~~~~~~~
.. autoclass:: WireByConcat

PlaneBySection
~~~~~~~~~~~~~~
.. autoclass:: PlaneBySection

PointsAlongShape
~~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongShape

PointsAlongEdge
~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongEdge

PointsAlongWire
~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongWire

FaceByPlane
~~~~~~~~~~~
.. autoclass:: FaceByPlane

WireByOffset
~~~~~~~~~~~~
.. autoclass:: WireByOffset

ShellByPipe
~~~~~~~~~~~
.. autoclass:: ShellByPipe

HalfspaceByShape
~~~~~~~~~~~~~~~~
.. autoclass:: HalfspaceByShape

WiresFromShape
~~~~~~~~~~~~~~
.. autoclass:: WiresFromShape

ShapeByOffset
~~~~~~~~~~~~~
.. autoclass:: ShapeByOffset

SolidsFromShapes
~~~~~~~~~~~~~~~~
.. autoclass:: SolidsFromShapes

Explore
-------
.. py:currentmodule:: afem.occ.explore

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
.. py:currentmodule:: afem.occ.modify

UnifyShape
~~~~~~~~~~
.. autoclass:: UnifyShape

FixShape
~~~~~~~~
.. autoclass:: FixShape

DivideClosedShape
~~~~~~~~~~~~~~~~~
.. autoclass:: DivideClosedShape

DivideC0Shape
~~~~~~~~~~~~~~~~~
.. autoclass:: DivideC0Shape

SplitWire
~~~~~~~~~
.. autoclass:: SplitWire

Boolean
-------
.. py:currentmodule:: afem.occ.bop

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

Properties
----------
Tools for shape properties (length, area, volume, etc.) will be here.

Check
-----

CheckShape
~~~~~~~~~~
.. autoclass:: afem.occ.check.CheckShape

Utilities
---------
.. automodule:: afem.occ.utils