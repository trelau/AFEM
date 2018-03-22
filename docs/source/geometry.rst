Geometry
========
The ``geometry`` package provides entities and tools for the creation and use of
what is commonly referred to as "construction geometry" or "reference geometry"
for both 2-D and 3-D domains. This package primarily wraps a number of
OpenCASCADE native types and tools in order to provide a more "Pythonic" user
interface. The entities and tools in the ``Geometry`` package do not cover
every OpenCASCADE type, but rather those frequently encountered during regular
AFEM use. The entities and tools can be imported by::

    from afem.geometry import *

Entities
--------
.. py:currentmodule:: afem.geometry.entities

Geometry2D
~~~~~~~~~~
.. autoclass:: Geometry2D

Point2D
~~~~~~~
.. autoclass:: Point2D

NurbsCurve2D
~~~~~~~~~~~~
.. autoclass:: NurbsCurve2D

Geometry
~~~~~~~~
.. autoclass:: Geometry

Point
~~~~~
.. autoclass:: Point

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

Circle
~~~~~~
.. autoclass:: Circle

Ellipse
~~~~~~~
.. autoclass:: Ellipse

NurbsCurve
~~~~~~~~~~
.. autoclass:: NurbsCurve

TrimmedCurve
~~~~~~~~~~~~
.. autoclass:: TrimmedCurve

Surface
~~~~~~~
.. autoclass:: Surface

Plane
~~~~~
.. autoclass:: Plane

NurbsSurface
~~~~~~~~~~~~
.. autoclass:: NurbsSurface

Create
------
.. py:currentmodule:: afem.geometry.create

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

CircleByNormal
~~~~~~~~~~~~~~
.. autoclass:: CircleByNormal

CircleByPlane
~~~~~~~~~~~~~
.. autoclass:: CircleByPlane

CircleBy3Points
~~~~~~~~~~~~~~~
.. autoclass:: CircleBy3Points

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

TrimmedCurveByParameters
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: TrimmedCurveByParameters

TrimmedCurveByCurve
~~~~~~~~~~~~~~~~~~~
.. autoclass:: TrimmedCurveByCurve

TrimmedCurveByPoints
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: TrimmedCurveByPoints

PlaneByNormal
~~~~~~~~~~~~~
.. autoclass:: PlaneByNormal

PlaneByAxes
~~~~~~~~~~~
.. autoclass:: PlaneByAxes

PlaneByPoints
~~~~~~~~~~~~~
.. autoclass:: PlaneByPoints

PlaneByApprox
~~~~~~~~~~~~~
.. autoclass:: PlaneByApprox

PlaneFromParameter
~~~~~~~~~~~~~~~~~~
.. autoclass:: PlaneFromParameter

PlaneByOrientation
~~~~~~~~~~~~~~~~~~
.. autoclass:: PlaneByOrientation

PlaneByCurveAndSurface
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlaneByCurveAndSurface

PlanesAlongCurveByNumber
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlanesAlongCurveByNumber

PlanesAlongCurveByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlanesAlongCurveByDistance

PlanesBetweenPlanesByNumber
~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlanesBetweenPlanesByNumber

PlanesBetweenPlanesByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlanesBetweenPlanesByDistance

PlanesAlongCurveAndSurfaceByDistance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: PlanesAlongCurveAndSurfaceByDistance

NurbsSurfaceByData
~~~~~~~~~~~~~~~~~~
.. autoclass:: NurbsSurfaceByData

NurbsSurfaceByInterp
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: NurbsSurfaceByInterp

NurbsSurfaceByApprox
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: NurbsSurfaceByApprox

Project
-------
.. py:currentmodule:: afem.geometry.project

PointProjector
~~~~~~~~~~~~~~
.. autoclass:: PointProjector

ProjectPointToCurve
~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectPointToCurve

ProjectPointToSurface
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectPointToSurface

CurveProjector
~~~~~~~~~~~~~~
.. autoclass:: CurveProjector

ProjectCurveToPlane
~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectCurveToPlane

ProjectCurveToSurface
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: ProjectCurveToSurface

Intersect
---------
.. py:currentmodule:: afem.geometry.intersect

CurveIntersector
~~~~~~~~~~~~~~~~
.. autoclass:: CurveIntersector

IntersectCurveCurve
~~~~~~~~~~~~~~~~~~~
.. autoclass:: IntersectCurveCurve

IntersectCurveSurface
~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: IntersectCurveSurface

SurfaceIntersector
~~~~~~~~~~~~~~~~~~
.. autoclass:: SurfaceIntersector

IntersectSurfaceSurface
~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: IntersectSurfaceSurface

Distance
--------
.. py:currentmodule:: afem.geometry.distance

DistancePointToCurve
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistancePointToCurve

DistancePointToSurface
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistancePointToSurface

DistanceCurveToCurve
~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistanceCurveToCurve

DistanceCurveToSurface
~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistanceCurveToSurface

DistanceSurfaceToSurface
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: DistanceSurfaceToSurface

Check
-----

CheckGeom
~~~~~~~~~
.. py:currentmodule:: afem.geometry.check
.. autoclass:: CheckGeom

Utilities
---------
.. automodule:: afem.geometry.utils