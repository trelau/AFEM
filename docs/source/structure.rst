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

Create
------

Create Part Module
~~~~~~~~~~~~~~~~~~
.. automodule:: afem.structure.creator

.. py:currentmodule:: afem.structure.create

CurvePartByShape
~~~~~~~~~~~~~~~~
.. autoclass:: CurvePartByShape

SurfacePartByShape
~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfacePartByShape

SparByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: SparByParameters

SparByPoints
~~~~~~~~~~~~
.. autoclass:: SparByPoints

SparByShape
~~~~~~~~~~~
.. autoclass:: SparByShape

SparBetweenGeom
~~~~~~~~~~~~~~~
.. autoclass:: SparBetweenGeom

SparsBetweenPlanes
~~~~~~~~~~~~~~~~~~
.. autoclass:: SparsBetweenPlanes

SparsAtShapes
~~~~~~~~~~~~~
.. autoclass:: SparsAtShapes

SparsAlongCurve
~~~~~~~~~~~~~~~
.. autoclass:: SparsAlongCurve

RibByParameters
~~~~~~~~~~~~~~~~
.. autoclass:: RibByParameters

RibByPoints
~~~~~~~~~~~~
.. autoclass:: RibByPoints

RibByShape
~~~~~~~~~~~
.. autoclass:: RibByShape

RibBetweenGeom
~~~~~~~~~~~~~~~
.. autoclass:: RibBetweenGeom

RibsBetweenPlanes
~~~~~~~~~~~~~~~~~~
.. autoclass:: RibsBetweenPlanes

RibsAtShapes
~~~~~~~~~~~~~
.. autoclass:: RibsAtShapes

RibsAlongCurve
~~~~~~~~~~~~~~~
.. autoclass:: RibsAlongCurve

BulkheadByShape
~~~~~~~~~~~~~~~
.. autoclass:: BulkheadByShape

FloorByShape
~~~~~~~~~~~~
.. autoclass:: FloorByShape

FrameByShape
~~~~~~~~~~~~
.. autoclass:: FrameByShape

FramesBetweenPlanes
~~~~~~~~~~~~~~~~~~~
.. autoclass:: FramesBetweenPlanes

FramesAtShapes
~~~~~~~~~~~~~~
.. autoclass:: FramesAtShapes

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

Modify
------

Tools
-----

.. autoclass:: afem.structure.tools.PartTools

Check
-----
.. autoclass:: afem.structure.checker.CheckPart