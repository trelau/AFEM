# AFEM — Airframe Finite Element Modeler
The AFEM application is an airframe modeling and FEM generation tool developed
by Laughlin Research, LLC. AFEM enables the use of high-order structural
analysis in the early phases of aircraft design.

Although AFEM targets airframe design and analysis applications, really only
the ``afem.structure`` and ``afem.oml`` packages directly support that cause.
All the other packages provide a more general and "Pythonic" set of entities
and tools that could potentially be used to develop applications in other
disciplines and/or domains.

## Technology Stack
The AFEM core technology stack includes:

* [Python](https://www.python.org/): The Python programming language enables
  rapid development and integration with other systems.

* [OpenCASCADE](https://www.opencascade.com): This mature library provides
  advanced geometric modeling and CAD functionality and is under active
  development.

* [Netgen](https://sourceforge.net/projects/netgen-mesher): This library
  enables advanced meshing capabilities including 3-D tetrahedral and 2-D
  unstructured quad-dominated surface meshing.

* [Salome Platform](http://www.salome-platform.org): The core meshing library
  from this open source application serves as the central component for
  AFEM's mesh generation capabilities.
  
* [pyOCCT](https://github.com/LaughlinResearch/pyOCCT): This open source
  project provides Python bindings to the OpenCASCADE and Salome Platform
  meshing libraries.

# Installation
AFEM is currently only supported for Windows 64-bit Python 3.5 and 3.6.
Cross-platform support is dependent upon prerequisites such as pyOCCT.
[Anaconda Python](https://www.anaconda.com/download/) or
[Miniconda](https://conda.io/miniconda.html) is recommended for installation
and regular use since many of the prerequisites are available via the Anaconda
Cloud.

It is recommended that a designated environment be created and used for AFEM.
An example of creating this environment for Anaconda Python within an Anaconda
command prompt is:

    conda create -n afem python=3.5

This will create an environment named "afem" with Python 3.5. Make sure this
environment is active when using AFEM. For Anaconda Python, activating this
environment may look like:

    activate afem

within an Anaconda command prompt.
 
The [pyOCCT](https://github.com/LaughlinResearch/pyOCCT) package developed by
Laughlin Research should now be installed. For supported platforms, installing
pyOCCT can be done by:

    conda install -c trelau pyocct

Other dependencies such as NumPy and SciPy can be installed as needed using
the conda package manager:

    conda install numpy scipy
    
or pip:

    pip install numpy scipy
    
A minimal graphical user interface requires the wxPython package which can be
installed by:

    conda install -c conda-forge wxpython

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
index.html file with a web browser.

# Getting Started
The best way to get started is to examine and run the files in the *examples/*
folder. Running the script:

    python structure_wingbox.py
    
should generate an image similar to the one shown below. Remember to make sure
the appropriate environment is active when using AFEM is applicable.

![wingbox](./docs/source/resources/wingbox.png)
