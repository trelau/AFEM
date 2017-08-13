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

Assembly
~~~~~~~~
.. autoclass:: afem.structure.assembly.Assembly

AssemblyAPI
~~~~~~~~~~~
.. autoclass:: afem.structure.assembly.AssemblyAPI

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

SparBySurface
~~~~~~~~~~~~~
.. autoclass:: SparBySurface

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

SparsAtShapes
~~~~~~~~~~~~~
.. autoclass:: SparsAtShapes

RibByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: RibByParameters

RibByPoints
~~~~~~~~~~~~
.. autoclass:: RibByPoints

RibBySurface
~~~~~~~~~~~~
.. autoclass:: RibBySurface

RibBetweenShapes
~~~~~~~~~~~~~~~~
.. autoclass:: RibBetweenShapes

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

RibsAtShapes
~~~~~~~~~~~~~
.. autoclass:: RibsAtShapes

BulkheadBySurface
~~~~~~~~~~~~~~~~~
.. autoclass:: BulkheadBySurface

FloorBySurface
~~~~~~~~~~~~~~
.. autoclass:: FloorBySurface

FrameByPlane
~~~~~~~~~~~~~
.. autoclass:: FrameByPlane

FramesBetweenPlanesByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: FramesBetweenPlanesByNumber

FramesBetweenPlanesByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: FramesBetweenPlanesByDistance

FramesByPlanes
~~~~~~~~~~~~~~
.. autoclass:: FramesByPlanes

SkinBySolid
~~~~~~~~~~~
.. autoclass:: SkinBySolid

SkinByBody
~~~~~~~~~~
.. autoclass:: SkinByBody

StringerBySpine
~~~~~~~~~~~~~~~
.. autoclass:: StringerBySpine

StringerBySection
~~~~~~~~~~~~~~~~~
.. autoclass:: StringerBySection

StringersBySections
~~~~~~~~~~~~~~~~~~~
.. autoclass:: StringersBySections

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

Modify
------

Check
-----
.. autoclass:: afem.structure.check.CheckPart

Utilities
---------
.. automodule:: afem.structure.utils
