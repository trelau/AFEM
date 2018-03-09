Exchange
========
The ``Exchange`` package provides a number of common tools for data exchange
including BREP (OpenCASCADE native format), STEP, IGES, and STL. The tools can
be imported by::

    from afem.exchange import *

Some specific use-case tools are provided including streamlined processing of
STEP files from NASA's OpenVSP program. These data exchange tools can be used
to import or export data to other CAD or CAE tools.

BREP
----
The ``afem.exchange.brep`` module contains two simple methods for reading and
writing OpenCASCADE BREP files.

.. automodule:: afem.exchange.brep

STEP
----
.. automodule:: afem.exchange.step

IGES
----
.. automodule:: afem.exchange.iges

STL
---
.. automodule:: afem.exchange.stl

OpenVSP
-------
.. automodule:: afem.exchange.vsp

NASTRAN
-------
.. automodule:: afem.exchange.nastran