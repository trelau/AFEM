ASAP
====
The Aerospace Structural Analysis Program (ASAP) is an airframe modeling and
FEM generation tool developed by Laughlin Research, LLC. ASAP enables the use
of high-order structural analysis in the early phases of aircraft design.

Installation
============
ASAP is currently developed for Python 3.5. Anaconda Python is recommended
for package management and since pre-built binaries are available for the
ASAP prerequisites using the Anaconda cloud (only for Windows 64-bit).


Prerequisites
-------------
ASAP relies on a number of LGPL licensed open-source tools, including:

    - `OpenCASCADE Community Edition <https://github.com/tpaviot/oce/releases/tag/OCE-0.17.2>`_

    - `pythonocc-core<https://github.com/trelau/pythonocc-core/tree/review/smesh-support>`_

    - `Netgen<https://github.com/trelau/netgen/tree/netgen4smesh>`_

    - `SMESH<https://github.com/trelau/smesh/tree/review/fc-smesh-771>`_

Pre-built binaries for these tools are available through the Anaconda cloud
for Python 3.5 Windows 64-bit. It is recommended that a designated environment
be created and used for ASAP. An example of creating this environment for
Anaconda Python within an Anaconda command prompt is:

    conda create -n asap python=3.5

This will create an environment named "asap" with Python 3.5. Make sure this
environment is active when using ASAP. For Anaconda Python, activating this
environment may look like:

    activate asap

within an Anaconda command prompt. At this point the prerequisites can be
installed using specified channels on the Anaconda cloud:

    conda install -c trelau -c oce -c dlr-sc pythonocc-core=smesh

This should automatically resolve all dependencies and install all the
required packages.

Installing ASAP
---------------
ASAP is a pure Python package and can be installed using the command:

    python setup.py develop

within the ASAP root folder. The "develop" option links to the source code
at runtime so changes in the source are reflected in any programs using ASAP.

Getting Started
---------------
The best way to get started is to examine and run the files in the examples and
test folders.