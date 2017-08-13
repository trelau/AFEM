Installation
============
This section explains how to install AFEM.

Supported Platforms
-------------------
Currently, only Python 3.5 Windows 64-bit is supported. Although AFEM is
Python 2/3 compatible, some core dependencies have thus far only been built
for Python 3.5 Windows 64-bit. Eventually more platforms will be supported.

Installing Python
-----------------
For new users, `installing Anaconda <https://www.continuum.io/downloads>`_ is
**highly recommended**. Anaconda conveniently installs Python and many other
commonly used packages for scientific computing. This installation
documentation assumes the use of Anaconda Python.

Many of the packages required by AFEM will be included in the default Anaconda
installation including NumPy and SciPy. See
`here <https://docs.continuum.io/anaconda/pkg-docs>`_ for a complete list
and more information about the Anaconda installation.

A specific environment for AFEM should be created using an Anaconda Command
Prompt::

    conda create -n afem python=3.5

The *-n* flag signifies the name and in this example it is *afem*, but it
can be anything. The *python=3.5* command tells Anaconda to create an
environment based on Python 3.5.

Make sure this environment is active when installing and using AFEM. For
Anaconda Python, activating this environment looks like::

    activate afem

within an Anaconda Command Prompt.

Installing Dependencies
-----------------------
AFEM relies on a number of LGPL licensed open-source tools, including:

    * `OpenCASCADE Community Edition <https://github.com/tpaviot/oce/releases/tag/OCE-0.18.1>`_
    * `pythonocc-core <https://github.com/trelau/pythonocc-core/releases/tag/0.18.2>`_
    * `Netgen <https://github.com/trelau/netgen/releases/tag/6.3>`_
    * `SMESH <https://github.com/trelau/smesh/releases/tag/7.7.2>`_

These dependencies can be installed using the Anaconda Cloud. From an Anaconda
Command Prompt with the desired environment activated::

    conda install -c trelau -c oce -c dlr-sc -c conda-forge pythonocc-core=0.18.2

This should automatically resolve any other dependencies.

AFEM also makes use of the NumPy and SciPy packages. Install these packages
to the designated environment by::

    conda install numpy scipy

At this point all dependencies should be installed.

Installing AFEM
---------------
AFEM is currently distributed in source form and can be installed using the
standard command::

    python setup.py install

This command should be performed within an Anaconda Command Prompt with the
designated environment activated.

After installation examples are available and unit tests can be performed::

    cd docs
    python run_tests.py

Failed tests can be reported to support@laughlinresearch.com