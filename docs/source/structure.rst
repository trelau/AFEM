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
.. inheritance-diagram:: CurvePart
   :parts: 1

.. autoclass:: CurvePart

Beam
~~~~
.. inheritance-diagram:: Beam
   :parts: 1

.. autoclass:: Beam

SurfacePart
~~~~~~~~~~~
.. inheritance-diagram:: SurfacePart
   :parts: 1

.. autoclass:: SurfacePart

WingPart
~~~~~~~~
.. inheritance-diagram:: WingPart
   :parts: 1

.. autoclass:: WingPart

Spar
~~~~
.. inheritance-diagram:: Spar
   :parts: 1

.. autoclass:: Spar

Rib
~~~
.. inheritance-diagram:: Rib
   :parts: 1

.. autoclass:: Rib

FuselagePart
~~~~~~~~~~~~
.. inheritance-diagram:: FuselagePart
   :parts: 1

.. autoclass:: FuselagePart

Bulkhead
~~~~~~~~
.. inheritance-diagram:: Bulkhead
   :parts: 1

.. autoclass:: Bulkhead

Floor
~~~~~
.. inheritance-diagram:: Floor
   :parts: 1

.. autoclass:: Floor

Frame
~~~~~
.. inheritance-diagram:: Frame
   :parts: 1

.. autoclass:: Frame

Skin
~~~~
.. inheritance-diagram:: Skin
   :parts: 1

.. autoclass:: Skin

Stiffener1D
~~~~~~~~~~~
.. inheritance-diagram:: Stiffener1D
   :parts: 1

.. autoclass:: Stiffener1D

Stiffener2D
~~~~~~~~~~~
.. inheritance-diagram:: Stiffener2D
   :parts: 1

.. autoclass:: Stiffener2D

Stringer
~~~~~~~~
.. inheritance-diagram:: Stringer
   :parts: 1

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

CutParts
~~~~~~~~
.. autoclass:: CutParts

SewSurfaceParts
~~~~~~~~~~~~~~~
.. autoclass:: SewSurfaceParts

SplitParts
~~~~~~~~~~
.. autoclass:: SplitParts

FuseAssemblies
~~~~~~~~~~~~~~
.. autoclass:: FuseAssemblies

Modify
------
.. py:currentmodule:: afem.structure.modify

DiscardByCref
~~~~~~~~~~~~~
.. autoclass:: DiscardByCref

Check
-----
.. autoclass:: afem.structure.check.CheckPart

Utilities
---------
.. automodule:: afem.structure.utils
