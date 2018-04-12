Structure
=========
**NOT YET COMPLETE***
This section describes the structure package. The entities and tools can be
imported by::

    from afem.structure import *

Entities
--------
.. py:currentmodule:: afem.structure.entities

Part
~~~~
.. autoclass:: Part

CurvePart
~~~~~~~~~
.. autoclass:: CurvePart

Beam1D
~~~~~~
.. autoclass:: Beam1D

SurfacePart
~~~~~~~~~~~
.. autoclass:: SurfacePart

WingPart
~~~~~~~~
.. autoclass:: WingPart

Spar
~~~~
.. autoclass:: Spar

Rib
~~~
.. autoclass:: Rib

FuselagePart
~~~~~~~~~~~~
.. autoclass:: FuselagePart

Bulkhead
~~~~~~~~
.. autoclass:: Bulkhead

Floor
~~~~~
.. autoclass:: Floor

Frame
~~~~~
.. autoclass:: Frame

Skin
~~~~
.. autoclass:: Skin

Stiffener1D
~~~~~~~~~~~
.. autoclass:: Stiffener1D

Stiffener2D
~~~~~~~~~~~
.. autoclass:: Stiffener2D

Stringer
~~~~~~~~
.. autoclass:: Stringer

Beam2D
~~~~~~
.. autoclass:: Beam2D

Group
~~~~~
.. autoclass:: afem.structure.group.Group

GroupAPI
~~~~~~~~
.. autoclass:: afem.structure.group.GroupAPI

Create
------
.. py:currentmodule:: afem.structure.create

PartBuilder
~~~~~~~~~~~
.. autoclass:: PartBuilder

PartsBuilder
~~~~~~~~~~~~
.. autoclass:: PartsBuilder

CurvePartByShape
~~~~~~~~~~~~~~~~
.. autoclass:: CurvePartByShape

Beam1DByShape
~~~~~~~~~~~~~
.. autoclass:: Beam1DByShape

Beam1DByCurve
~~~~~~~~~~~~~
.. autoclass:: Beam1DByCurve

Beam1DByPoints
~~~~~~~~~~~~~~
.. autoclass:: Beam1DByPoints

SurfacePartByShape
~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartByShape

SurfacePartByParameters
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartByParameters

SurfacePartByPoints
~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartByPoints

SurfacePartByEnds
~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartByEnds

SurfacePartBetweenShapes
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartBetweenShapes

SurfacePartsBetweenPlanesByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartsBetweenPlanesByNumber

SurfacePartsBetweenPlanesByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartsBetweenPlanesByDistance

SurfacePartsAlongCurveByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartsAlongCurveByNumber

SurfacePartsAlongCurveByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartsAlongCurveByDistance

SparByShape
~~~~~~~~~~~
.. autoclass:: SparByShape

SparByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: SparByParameters

SparByPoints
~~~~~~~~~~~~
.. autoclass:: SparByPoints

SparByEnds
~~~~~~~~~~
.. autoclass:: SparByEnds

SparBetweenShapes
~~~~~~~~~~~~~~~~~
.. autoclass:: SparBetweenShapes

SparsBetweenPlanesByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SparsBetweenPlanesByNumber

SparsBetweenPlanesByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SparsBetweenPlanesByDistance

SparsAlongCurveByNumber
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SparsAlongCurveByNumber

SparsAlongCurveByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: SparsAlongCurveByDistance

RibByShape
~~~~~~~~~~
.. autoclass:: RibByShape

RibByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: RibByParameters

RibByPoints
~~~~~~~~~~~~
.. autoclass:: RibByPoints

RibBetweenShapes
~~~~~~~~~~~~~~~~
.. autoclass:: RibBetweenShapes

RibByOrientation
~~~~~~~~~~~~~~~~
.. autoclass:: RibByOrientation

RibsBetweenPlanesByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RibsBetweenPlanesByNumber

RibsBetweenPlanesByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RibsBetweenPlanesByDistance

RibsAlongCurveByNumber
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RibsAlongCurveByNumber

RibsAlongCurveByDistance
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RibsAlongCurveByDistance

RibsAlongCurveAndSurfaceByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: RibsAlongCurveAndSurfaceByDistance

BulkheadByShape
~~~~~~~~~~~~~~~
.. autoclass:: BulkheadByShape

FloorByShape
~~~~~~~~~~~~
.. autoclass:: FloorByShape

FrameByPlane
~~~~~~~~~~~~~
.. autoclass:: FrameByPlane

FramesByPlanes
~~~~~~~~~~~~~~
.. autoclass:: FramesByPlanes

FramesBetweenPlanesByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: FramesBetweenPlanesByNumber

FramesBetweenPlanesByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: FramesBetweenPlanesByDistance

SkinBySolid
~~~~~~~~~~~
.. autoclass:: SkinBySolid

SkinByBody
~~~~~~~~~~
.. autoclass:: SkinByBody

StringerByShape
~~~~~~~~~~~~~~~
.. autoclass:: StringerByShape

Beam2DBySweep
~~~~~~~~~~~~~
.. autoclass:: Beam2DBySweep

Explore
-------

Join
----
.. py:currentmodule:: afem.structure.join

FuseSurfaceParts
~~~~~~~~~~~~~~~~
.. autoclass:: FuseSurfaceParts

FuseSurfacePartsByCref
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: FuseSurfacePartsByCref

CutParts
~~~~~~~~
.. autoclass:: CutParts

SewSurfaceParts
~~~~~~~~~~~~~~~
.. autoclass:: SewSurfaceParts

SplitParts
~~~~~~~~~~
.. autoclass:: SplitParts

FuseGroups
~~~~~~~~~~
.. autoclass:: FuseGroups

Modify
------
.. py:currentmodule:: afem.structure.modify

DiscardByCref
~~~~~~~~~~~~~
.. autoclass:: DiscardByCref

Fix
---
.. py:currentmodule:: afem.structure.fix

FixGroup
~~~~~~~~
.. autoclass:: FixGroup

Check
-----
.. autoclass:: afem.structure.check.CheckPart

Utilities
---------
.. automodule:: afem.structure.utils
