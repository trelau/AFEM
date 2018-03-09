# AFEM — Airframe Finite Element Modeler
The AFEM application is an airframe modeling and FEM generation tool developed
by Laughlin Research, LLC. AFEM enables the use of high-order structural
analysis in the early phases of aircraft design.

## Technology Stack
The AFEM core technology stack includes:

* [OpenCASCADE](https://www.opencascade.com): This mature library provides
  advanced geometric modeling and CAD functionality and is under active
  development.

* [Netgen](https://sourceforge.net/projects/netgen-mesher): This library
  enables advanced meshing capabilities including 3-D tetrahedral and 2-D
  unstructured quad-dominated surface meshing.

* [Salome Platform](http://www.salome-platform.org): The core meshing library
  from this open source application serves as the central component for
  AFEM's mesh generation capabilities.
  
* [pyOCCT](https://github.com/LaughlinResearch/pyOCCT): This open source project
  provides Python bindings to the OpenCASCADE and Salome Platform meshing
  libraries.

# Installation
AFEM is currently only supported for Windows 64-bit Python 3.5. [Anaconda
Python](https://www.anaconda.com/download/) is recommended for package
management and since prebuilt binaries are available for the prerequisites.

It is recommended that a designated environment be created and used for AFEM. An
example of creating this environment for Anaconda Python within an Anaconda
command prompt is:

    conda create -n afem python=3.5

This will create an environment named "afem" with Python 3.5. Make sure this
environment is active when using AFEM. For Anaconda Python, activating this
environment may look like:

    activate afem

within an Anaconda command prompt.
 
The [pyOCCT](https://github.com/LaughlinResearch/pyOCCT) package developed by
Laughlin Research should now be installed. A Python wheel is available for
Windows 64-bit Python 3.5 that contains the pyOCCT package as well as all
necessary dependencies and can be installed using
[pip](https://pypi.python.org/pypi/pip/). From the directory where the pyOCCT
wheel is located the installation command may look like:

    pip install OCCT-0.0.1-cp35-none-win_amd64.whl
    
Note that the wheel filename may be different depending on version and platform.

Other dependencies such as NumPy and SciPy can be installed as needed using
the conda package manager:

    conda install numpy scipy
    
or pip:

    pip install numpy scipy
    
A minimal graphical user interface requires the PySide package which can be
installed by:

    conda install -c conda-forge pyside=1.2.4

# Installing AFEM
Be sure to activate the designed AFEM environment before installation.

AFEM is a pure Python package and can be installed using the command:

    python setup.py develop

within the AFEM root folder. The ``develop`` option links to the source code
at runtime so changes in the source are reflected in any programs using AFEM.
The regular installation command:

    python setup.py install
    
can be used to actually install the AFEM package into the Python site-packages
directory.

Installation files can be cleaned by:

    conda clean -a

# Building Documentation
The documentation can be built from sources using sphinx. Install sphinx and
sphinx_rtd_theme in the desired conda environment by:

    conda install sphinx sphinx_rtd_theme
    
Then navigate to the *docs/* folder and run:

    make html

This should build html documentation in a *docs/build/html* folder. Open the 
afem.html file with a web browser.

# Getting Started
The best way to get started is to examine and run the files in the examples and
test folders.

# Notice
Copyright (c) 2018, Laughlin Research, LLC

Terms of Use:

The AFEM Code, including its source code and related software
documentation (collectively, the "AFEM Code"), as distributed herein
and as may be subsequently revised, in whole and in part, is for
government use only pursuant to development agreements between NASA,
Georgia Institute of Technology, and Laughlin Research, LLC. At the
time of distribution hereof, none of the AFEM Code is believed or
intended to be open source. Disclosure of the AFEM Code is strictly
subject to one or more restrictive covenants, including
non-disclosure and non-circumvention covenants, and any use of the
whole or a part of the AFEM Code constitutes acknowledgement and
acceptance of said covenants. Any unauthorized use, disclosure,
and/or sale of the AFEM Code or any portion thereof may be actionable
under current law.
