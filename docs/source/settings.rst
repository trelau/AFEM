Settings
========
Default settings for the AFEM program are provided in the :class:`.Settings`
class and accessed as class variables. These settings can be adjusted during
use and are passed to many of the classes and methods. The :class:`.Settings`
class can be accessed by::

    from afem.config import Settings

A logging utility is used to provide useful information during program
execution. A file with the name *afem.log* will be automatically created
wherever the main script is executed and whose contents will be dependent on
the logging level. In order to output the logging content to the command window
the following method should be called before the main script begins::

    Settings.log_to_console()

The logging content will now be displayed in the command window as well as
output to the log file.

Perhaps the setting with the most implication is what units are set for
OpenCASCADE. This is especially critical during data exchange activities
(e.g., STEP translation). The units can be set using the class method::

    Settings.set_units('in')

Inches ('in'), feet ('ft'), meters ('m'), and millimeters ('mm') are available.
OpenCASCADE uses millimeters by default, but AFEM should use inches as its
default setting when units are relevant.

.. autoclass:: afem.config.Settings
