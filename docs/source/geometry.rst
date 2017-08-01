Geometry
========
This section describes the geometry package. The entities and tools can be 
imported by::

    from afem.geometry import *

Entities
--------
.. py:currentmodule:: afem.geometry.entities

Geometry
~~~~~~~~
.. autoclass:: Geometry

Point
~~~~~
.. autoclass:: Point

Point2D
~~~~~~~
.. autoclass:: Point2D

Direction
~~~~~~~~~
.. autoclass:: Direction

Vector
~~~~~~
.. autoclass:: Vector

Axis1
~~~~~
.. autoclass:: Axis1

Axis3
~~~~~
.. autoclass:: Axis3

Curve
~~~~~
.. autoclass:: Curve

Line
~~~~
.. autoclass:: Line

NurbsCurve
~~~~~~~~~~
.. autoclass:: NurbsCurve

NurbsCurve2D
~~~~~~~~~~~~
.. autoclass:: NurbsCurve2D

Surface
~~~~~~~
.. autoclass:: Surface

Plane
~~~~~
.. autoclass:: Plane

NurbsSurface
~~~~~~~~~~~~
.. autoclass:: NurbsSurface

Bounding Box
~~~~~~~~~~~~
.. autoclass:: BBox

Create
------
.. py:currentmodule:: afem.geometry.create

CreateGeom
~~~~~~~~~~
.. autoclass:: CreateGeom

CreatedPoints
~~~~~~~~~~~~~
.. autoclass:: CreatedPoints

GeomBuilder
~~~~~~~~~~~
.. autoclass:: GeomBuilder

PointByXYZ
~~~~~~~~~~
.. autoclass:: PointByXYZ

PointByArray
~~~~~~~~~~~~
.. autoclass:: PointByArray

PointFromParameter
~~~~~~~~~~~~~~~~~~
.. autoclass:: PointFromParameter

PointsAlongCurveByNumber
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongCurveByNumber

PointsAlongCurveByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongCurveByDistance

DirectionByXYZ
~~~~~~~~~~~~~~
.. autoclass:: DirectionByXYZ

DirectionByArray
~~~~~~~~~~~~~~~~
.. autoclass:: DirectionByArray

DirectionByPoints
~~~~~~~~~~~~~~~~~
.. autoclass:: DirectionByPoints

VectorByXYZ
~~~~~~~~~~~
.. autoclass:: VectorByXYZ

VectorByArray
~~~~~~~~~~~~~
.. autoclass:: VectorByArray

VectorByPoints
~~~~~~~~~~~~~~
.. autoclass:: VectorByPoints

LineByVector
~~~~~~~~~~~~
.. autoclass:: LineByVector

LineByPoints
~~~~~~~~~~~~
.. autoclass:: LineByPoints

NurbsCurveByData
~~~~~~~~~~~~~~~~
.. autoclass:: NurbsCurveByData

NurbsCurveByInterp
~~~~~~~~~~~~~~~~~~
.. autoclass:: NurbsCurveByInterp

NurbsCurveByApprox
~~~~~~~~~~~~~~~~~~
.. autoclass:: NurbsCurveByApprox

NurbsCurveByPoints
~~~~~~~~~~~~~~~~~~
.. autoclass:: NurbsCurveByPoints

CurveByUIso
~~~~~~~~~~~
.. autoclass:: CurveByUIso

CurveByVIso
~~~~~~~~~~~
.. autoclass:: CurveByVIso

PlaneByNormal
~~~~~~~~~~~~~
.. autoclass:: PlaneByNormal

PlaneByAxes
~~~~~~~~~~~
.. autoclass:: PlaneByAxes

PlaneByPoints
~~~~~~~~~~~~~
.. autoclass:: PlaneByPoints

PlaneByFit
~~~~~~~~~~
.. autoclass:: PlaneByFit

PlanesAlongCurve
~~~~~~~~~~~~~~~~
.. autoclass:: PlanesAlongCurve

PlanesBetweenPlanes
~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlanesBetweenPlanes

SurfaceByData
~~~~~~~~~~~~~
.. autoclass:: SurfaceByData

SurfaceByInterp
~~~~~~~~~~~~~~~
.. autoclass:: SurfaceByInterp

SurfaceByFit
~~~~~~~~~~~~
.. autoclass:: SurfaceByFit

Project
-------
.. py:currentmodule:: afem.geometry.project

ProjectGeom
~~~~~~~~~~~
.. autoclass:: ProjectGeom

ProjectPointToCurve
~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectPointToCurve
   :inherited-members:

ProjectPointToSurface
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectPointToSurface
   :inherited-members:

ProjectCurveToPlane
~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectCurveToPlane
   :inherited-members:

ProjectCurveToSurface
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectCurveToSurface
   :inherited-members:

Intersect
---------
.. py:currentmodule:: afem.geometry.intersect

IntersectGeom
~~~~~~~~~~~~~
.. autoclass:: IntersectGeom

IntersectCurveCurve
~~~~~~~~~~~~~~~~~~~
.. autoclass:: IntersectCurveCurve
   :inherited-members:

IntersectCurveSurface
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: IntersectCurveSurface
   :inherited-members:

IntersectSurfaceSurface
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: IntersectSurfaceSurface
   :inherited-members:

IntersectError
~~~~~~~~~~~~~~
.. autoclass:: IntersectError

Check
-----

CheckGeom
~~~~~~~~~
.. py:currentmodule:: afem.geometry.check
.. autoclass:: CheckGeom

Utilities
---------
.. automodule:: afem.geometry.utils