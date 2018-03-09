OML
===
The ``OML`` package provides a minimal set of entities and tools intended to
make defining structural components easier and more efficient. Many of the
tools used to create parts require at minimum a solid defining the common
material of the part given its basis shape. Some tools are more automated and
parametric and therefore more information, like a reference surface or curve, is
required. The ``OML`` package, and more specifically the ``Body`` class, are
intended define a single type that contains all necessary information. The
entities and tools can be imported by::

    from afem.oml import *

Entities
--------
.. py:currentmodule:: afem.oml.entities

Body
~~~~
.. autoclass:: Body

Check
-----

CheckOML
~~~~~~~~
.. autoclass:: afem.oml.check.CheckOML