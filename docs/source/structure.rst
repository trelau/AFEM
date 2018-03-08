Structure
=========
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

Beam
~~~~
.. autoclass:: Beam

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

Group
~~~~~
.. autoclass:: afem.structure.group.Group

GroupAPI
~~~~~~~~
.. autoclass:: afem.structure.group.GroupAPI

Create
------
.. py:currentmodule:: afem.structure.create

CurvePartByShape
~~~~~~~~~~~~~~~~
.. autoclass:: CurvePartByShape

BeamByShape
~~~~~~~~~~~
.. autoclass:: BeamByShape

BeamByCurve
~~~~~~~~~~~
.. autoclass:: BeamByCurve

BeamByPoints
~~~~~~~~~~~~
.. autoclass:: BeamByPoints

SurfacePartByShape
~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartByShape

SparByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: SparByParameters

SparByPoints
~~~~~~~~~~~~
.. autoclass:: SparByPoints

SparByEnds
~~~~~~~~~~
.. autoclass:: SparByEnds

SparBySurface
~~~~~~~~~~~~~
.. autoclass:: SparBySurface

SparByShape
~~~~~~~~~~~
.. autoclass:: SparByShape

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

RibByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: RibByParameters

RibByPoints
~~~~~~~~~~~~
.. autoclass:: RibByPoints

RibBySurface
~~~~~~~~~~~~
.. autoclass:: RibBySurface

RibByShape
~~~~~~~~~~
.. autoclass:: RibByShape

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

BulkheadBySurface
~~~~~~~~~~~~~~~~~
.. autoclass:: BulkheadBySurface

FloorBySurface
~~~~~~~~~~~~~~
.. autoclass:: FloorBySurface

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
