Settings
========
Default settings for the AFEM program are provided in the Settings class and
accessed as class variables. These settings can be adjusted during use and are
passed to many of the classes and methods. The Settings class can be
accessed by::

    from afem.config import Settings

A logging utility is used to provide useful information during program
execution. A file with the name *afem.log* will be automatically created
wherever the main script is executed and whose contents will be dependent on
the logging level. The Settings class provides a method to set the desired
output level of the logging utility.

.. autoclass:: afem.config.Settings
