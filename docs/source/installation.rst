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
The **pyOCCT** package developed by Laughlin Research should now be installed.
The source code, pre-built binaries for Windows Python 3.5, and installation
instructions can be found
[here](https://github.gatech.edu/LaughlinResearch/pyOCCT).

Other dependencies such as NumPy and SciPy can be installed as needed using
the conda package manager::

    conda install numpy scipy

Installing AFEM
---------------
Be sure to activate the designed AFEM environment before installation.

AFEM is a pure Python package and can be installed using the command:

    python setup.py develop

within the AFEM root folder. The ``develop`` option links to the source code
at runtime so changes in the source are reflected in any programs using AFEM.
The regular installation command

    python setup.py install

can be used to actually install the AFEM package into the Python root directory.

Installation files can be cleaned by:

    conda clean -a

Building Documentation
----------------------
The documentation can be built from sources using sphinx. Install sphinx and
sphinx_rtd_theme in the desired conda environment by::

    conda install sphinx sphinx_rtd_theme

Then navigate to the docs/ folder and run:

    make html

This should build html documentation in a docs/build/html folder. Open the
afem.html file with a web browser.
