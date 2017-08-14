FEM
===
This section describes the finite element modeling package. The entities and
tools can be imported by::

    from afem.fem import *

Entities
--------

Node
~~~~
.. autoclass:: afem.fem.nodes.Node

.. py:currentmodule:: afem.fem.elements

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
.. py:currentmodule:: afem.fem.hypotheses

Hypothesis
~~~~~~~~~~
.. autoclass:: Hypothesis

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

Meshes
------
.. py:currentmodule:: afem.fem.meshes

Mesh
~~~~
.. autoclass:: Mesh

SubMesh
~~~~~~~
.. autoclass:: SubMesh

MeshAPI
~~~~~~~
.. autoclass:: MeshAPI