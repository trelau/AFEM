Geometry
========
This section describes the geometry package. The entities and tools can be 
imported by::

    from afem.geometry import *

Entities
--------

Geometry
~~~~~~~~
.. autoclass:: afem.geometry.base.Geometry

Point
~~~~~
.. autoclass:: afem.geometry.points.Point

Point2D
~~~~~~~
.. autoclass:: afem.geometry.points.Point2D

Direction
~~~~~~~~~
.. autoclass:: afem.geometry.vectors.Direction

Vector
~~~~~~
.. autoclass:: afem.geometry.vectors.Vector

Axis1
~~~~~
.. autoclass:: afem.geometry.axes.Axis1

Axis3
~~~~~
.. autoclass:: afem.geometry.axes.Axis3

Curve
~~~~~
.. autoclass:: afem.geometry.curves.Curve

Line
~~~~
.. autoclass:: afem.geometry.curves.Line

NurbsCurve
~~~~~~~~~~
.. autoclass:: afem.geometry.curves.NurbsCurve

NurbsCurve2D
~~~~~~~~~~~~
.. autoclass:: afem.geometry.curves.NurbsCurve2D

Surface
~~~~~~~
.. autoclass:: afem.geometry.surfaces.Surface

Plane
~~~~~
.. autoclass:: afem.geometry.surfaces.Plane

NurbsSurface
~~~~~~~~~~~~
.. autoclass:: afem.geometry.surfaces.NurbsSurface

Bounding Box
~~~~~~~~~~~~
.. autoclass:: afem.geometry.bbox.BBox

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

PointFromOther
~~~~~~~~~~~~~~
.. autoclass:: PointFromOther

PointsAlongCurve
~~~~~~~~~~~~~~~~
.. autoclass:: PointsAlongCurve

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

CurveByData
~~~~~~~~~~~
.. autoclass:: CurveByData

CurveByInterp
~~~~~~~~~~~~~
.. autoclass:: CurveByInterp

CurveByFit
~~~~~~~~~~
.. autoclass:: CurveByFit

CurveByPoints
~~~~~~~~~~~~~
.. autoclass:: CurveByPoints

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