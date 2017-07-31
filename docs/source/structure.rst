Structure
=========
This section describes the structure package. The entities and tools can be
imported by::

    from afem.structure import *

Entities
--------

Part
~~~~
.. autoclass:: afem.structure.part.Part

CurvePart
~~~~~~~~~
.. autoclass:: afem.structure.curve_part.CurvePart

Beam
~~~~
.. autoclass:: afem.structure.beam.Beam

SurfacePart
~~~~~~~~~~~
.. autoclass:: afem.structure.surface_part.SurfacePart

WingPart
~~~~~~~~
.. autoclass:: afem.structure.wing_part.WingPart

Spar
~~~~
.. autoclass:: afem.structure.spar.Spar

Rib
~~~
.. autoclass:: afem.structure.rib.Rib

FuselagePart
~~~~~~~~~~~~
.. autoclass:: afem.structure.fuselage_part.FuselagePart

Bulkhead
~~~~~~~~
.. autoclass:: afem.structure.bulkhead.Bulkhead

Floor
~~~~~
.. autoclass:: afem.structure.floor.Floor

Frame
~~~~~
.. autoclass:: afem.structure.frame.Frame

Skin
~~~~
.. autoclass:: afem.structure.skin.Skin

Stiffener1D
~~~~~~~~~~~
.. autoclass:: afem.structure.stiffeners.Stiffener1D

Stiffener2D
~~~~~~~~~~~
.. autoclass:: afem.structure.stiffeners.Stiffener2D

Stringer
~~~~~~~~
.. autoclass:: afem.structure.stringer.Stringer

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