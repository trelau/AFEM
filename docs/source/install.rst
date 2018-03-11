Installation
============
This section explains how to install AFEM and required dependencies.

Currently, only Python 3.5 Windows 64-bit is supported. Although AFEM is
Python 2/3 compatible, some core dependencies have thus far only been built
for Python 3.5 Windows 64-bit. More platforms will be supported in the future.

Installing Python
-----------------
For new users, `installing Anaconda <https://www.continuum.io/downloads>`_ is
**highly recommended**. Anaconda conveniently installs Python and many other
commonly used packages for scientific computing. This installation
documentation assumes the use of Anaconda Python. Otherwise, the official
`Python <https://www.python.org/downloads/>`_ distribution can be used.

Many of the packages required by AFEM will be included in the default Anaconda
installation including NumPy and SciPy. See
`here <https://docs.continuum.io/anaconda/pkg-docs>`_ for a complete list
and more information about the Anaconda installation.

A specific environment for AFEM should be created using an Anaconda Command
Prompt::

    conda create -n afem python=3.5

The ``-n`` flag signifies the name and in this example it is *afem*, but it
can be anything. The ``python=3.5`` command tells Anaconda to create an
environment based on Python 3.5.

Make sure this environment is active when installing and using AFEM. For
Anaconda Python, activating this environment looks like::

    activate afem

within an Anaconda Command Prompt.

Installing Dependencies
-----------------------
The `pyOCCT <https://github.com/LaughlinResearch/pyOCCT>`_ package developed by
Laughlin Research should now be installed. A Python wheel is available for
Windows 64-bit Python 3.5 that contains the pyOCCT package as well as all
necessary dependencies and can be installed using
`pip <https://pypi.python.org/pypi/pip/>`_. From the directory where the pyOCCT
wheel is located the installation command may look like::

    pip install OCCT-0.0.1-cp35-none-win_amd64.whl

Note that the wheel filename may be different depending on version and platform.

Other dependencies such as NumPy and SciPy can be installed as needed using
the conda package manager::

    conda install numpy scipy

or pip::

    pip install numpy scipy

A minimal graphical user interface requires the wxPython package which can be
installed by::

    conda install -c conda-forge wxpython

or with pip::

    pip install wxpython

Installing AFEM
---------------
Be sure to activate the designed AFEM environment before installation. AFEM is a
pure Python package and can be installed using the command::

    python setup.py develop

within the AFEM root folder. The ``develop`` option links to the source code
at runtime so changes in the source are reflected in any programs using AFEM.
The regular installation command::

    python setup.py install

can be used to actually install the AFEM package into the Python *site-packages*
folder.

Installation files can be cleaned by::

    conda clean -a

If installing AFEM from a prebuilt wheel, use pip::

    pip install afem-0.2.0-cp35-none-win_amd64.whl

from the folder where the wheel is located. Be sure that the appropriate
environment is active before running this command.

Building Documentation
----------------------
The documentation can be built from sources using sphinx. Install sphinx and
sphinx_rtd_theme in the desired conda environment by::

    conda install sphinx sphinx_rtd_theme

Then navigate to the *docs/* folder and run::

    make html

This should build html documentation in a *docs/build/html* folder. Open the
afem.html file with a web browser.
