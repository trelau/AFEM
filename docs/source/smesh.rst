SMESH
=====
This section describes the meshing package. The entities and tools can be
imported by::

    from afem.smesh import *

This package leverages the standalone SMESH package extracted from the Salome
Platform.

Entities
--------
.. py:currentmodule:: afem.smesh.entities

Node
~~~~
.. autoclass:: Node

Element
~~~~~~~
.. autoclass:: Element

Hypotheses
----------
.. py:currentmodule:: afem.smesh.hypotheses

Hypothesis
~~~~~~~~~~
.. autoclass:: Hypothesis

Algorithm
~~~~~~~~~
.. autoclass:: Algorithm

Regular1D
~~~~~~~~~
.. autoclass:: Regular1D

CompositeSide1D
~~~~~~~~~~~~~~~
.. autoclass:: CompositeSide1D

MaxLength1D
~~~~~~~~~~~
.. autoclass:: MaxLength1D

LocalLength1D
~~~~~~~~~~~~~
.. autoclass:: LocalLength1D

NumberOfSegments1D
~~~~~~~~~~~~~~~~~~
.. autoclass:: NumberOfSegments1D

Adaptive1D
~~~~~~~~~~
.. autoclass:: Adaptive1D

Deflection1D
~~~~~~~~~~~~
.. autoclass:: Deflection1D

QuadrangleAlgo2D
~~~~~~~~~~~~~~~~
.. autoclass:: QuadrangleAlgo2D

QuadrangleHypo2D
~~~~~~~~~~~~~~~~
.. autoclass:: QuadrangleHypo2D

NetgenAlgo2D
~~~~~~~~~~~~
.. autoclass:: NetgenAlgo2D

NetgenAlgoOnly2D
~~~~~~~~~~~~~~~~
.. autoclass:: NetgenAlgoOnly2D

NetgenHypo2D
~~~~~~~~~~~~
.. autoclass:: NetgenHypo2D

NetgenSimple2D
~~~~~~~~~~~~~~
.. autoclass:: NetgenSimple2D

MeshGemsAlgo2D
~~~~~~~~~~~~~~
.. autoclass:: MeshGemsAlgo2D

MeshGemsHypo2D
~~~~~~~~~~~~~~
.. autoclass:: MeshGemsHypo2D

Meshes
------
.. py:currentmodule:: afem.smesh.meshes

MeshGen
~~~~~~~
.. autoclass:: MeshGen

Mesh
~~~~
.. autoclass:: Mesh

MeshDS
~~~~~~
.. autoclass:: MeshDS

SubMesh
~~~~~~~
.. autoclass:: SubMesh

SubMeshDS
~~~~~~~~~
.. autoclass:: SubMeshDS

Utilities
---------
.. py:currentmodule:: afem.smesh.utils

MeshEditor
~~~~~~~~~~
.. autoclass:: MeshEditor

MeshHelper
~~~~~~~~~~
.. autoclass:: MeshHelper
