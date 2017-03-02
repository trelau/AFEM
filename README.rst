ASAP
====
The Aerospace Structural Analysis Program (ASAP) is an airframe structural
design and modeling tool developed by Laughlin Research, LLC. The primary
goal of ASAP is to enable high-order structural modeling and analysis in the
early phases of aircraft design.

Installing
----------
ASAP can be installed using the command:

    python setup.py develop

ASAP relies on PythonOCC, SMESH, and Netgen which can be found here:

https://github.com/trelau

and installed using conda with:

    conda install -c trelau pythonocc-core

This relies on OpenCASCADE Community Edition 0.17.2 (OCE) which can be
installed with:

    conda install -c oce oce=0.17.2

Getting Started
---------------
The way to get started is to examine and run the files in the examples and
test folders.