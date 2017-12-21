Mesh
====
This section describes the meshing package. The entities and tools can be
imported by::

    from afem.mesh import *

Entities
--------

Node
~~~~
.. autoclass:: afem.mesh.nodes.Node

.. py:currentmodule:: afem.mesh.elements

Element
~~~~~~~
.. autoclass:: Element

Elm0D
~~~~~
.. autoclass:: Elm0D

Elm1D
~~~~~
.. autoclass:: Elm1D

Elm2D
~~~~~
.. autoclass:: Elm2D

Hypotheses
----------
.. py:currentmodule:: afem.mesh.hypotheses

Hypothesis
~~~~~~~~~~
.. autoclass:: Hypothesis

Algorithm
~~~~~~~~~
.. autoclass:: Algorithm

Regular1D
~~~~~~~~~
.. autoclass:: Regular1D

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

NetgenHypothesis
~~~~~~~~~~~~~~~~
.. autoclass:: NetgenHypothesis

NetgenSimple2D
~~~~~~~~~~~~~~
.. autoclass:: NetgenSimple2D

NetgenAlgo2D
~~~~~~~~~~~~
.. autoclass:: NetgenAlgo2D

NetgenAlgoOnly2D
~~~~~~~~~~~~~~~~
.. autoclass:: NetgenAlgoOnly2D

QuadrangleParams2D
~~~~~~~~~~~~~~~~~~
.. autoclass:: QuadrangleParams2D

Quadrangle2D
~~~~~~~~~~~~
.. autoclass:: Quadrangle2D

BlsurfAlgo
~~~~~~~~~~
.. autoclass:: BlsurfAlgo

BlsurfHypothesis
~~~~~~~~~~~~~~~~
.. autoclass:: BlsurfHypothesis

HypothesisAPI
~~~~~~~~~~~~~
.. autoclass:: HypothesisAPI

Meshes
------
.. py:currentmodule:: afem.mesh.meshes

Mesh
~~~~
.. autoclass:: Mesh

SubMesh
~~~~~~~
.. autoclass:: SubMesh

MeshAPI
~~~~~~~
.. autoclass:: MeshAPI