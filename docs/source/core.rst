Core
====
The ``afem.core`` package includes the ``ShapeHolder`` class which is a base
type meant to hold convenient properties and methods for the ``Body`` and
``Part`` types. A reference curve and surface can be set in this class as
needed by the user and a number of streamlined methods are available for
creating additional reference geometry. Included are methods to create points
and planes along the reference curve as well as projecting points to either the
reference curve or surface.

Entities
--------
.. py:currentmodule:: afem.core.entities

ShapeHolder
~~~~~~~~~~~
.. autoclass:: ShapeHolder
