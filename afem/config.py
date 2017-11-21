import logging
import sys

# Initialize logger.
with open('afem.log', 'w') as log:
    log.write('-----------------------------\n')
    log.write('AFEM LOGGING FILE INITIALIZED\n')
    log.write('-----------------------------\n')
logger = logging.getLogger('afem')
logger.setLevel(logging.INFO)
_fh = logging.FileHandler('afem.log')
_fmt = logging.Formatter('%(levelname)s: %(message)s')
_fh.setFormatter(_fmt)
logger.addHandler(_fh)

# Dictionary for units
units_dict = {'i': 'INCH',
              'in': 'INCH',
              'inch': 'INCH',
              'inches': 'INCH',
              'f': 'FOOT',
              'ft': 'FOOT',
              'foot': 'FOOT',
              'feet': 'FOOT',
              'm': 'M',
              'meter': 'M',
              'meters': 'M',
              'mm': 'MM',
              'millimeter': 'MM',
              'millimeters': 'MM'}

log_dict = {'debug': logging.DEBUG,
            'info': logging.INFO,
            'warning': logging.WARNING,
            'error': logging.ERROR,
            'critical': logging.CRITICAL}

__all__ = ["Settings"]


class Settings(object):
    """
    Settings.

    :var str units: The default units ('in', 'ft', 'm', 'mm'). The default
        value is inches.
    """
    # Class variables for settings
    units = 'INCH'

    @staticmethod
    def log_to_console():
        """
        Option to add a stream handler to the main file logger to log to
        stdout.

        :return: None.
        """
        chdlr = logging.StreamHandler(sys.stdout)
        logger.addHandler(chdlr)

    @staticmethod
    def set_loggging_level(level='info'):
        """
        Set the logging level.

        :param str level: Logging level ('debug', 'info', 'warning', 'error',
            'critical').

        :return: None.

        :raise KeyError: If the given level is not supported.

        .. note::
            Note all logging levels are currently utilized in the AFEM program.
            The default value of 'info' is recommended for general use.
        """
        level = level.lower()
        logger.setLevel(log_dict[level])

    @classmethod
    def set_units(cls, units='in'):
        """
        Set units.
        
        :param str units: The units ('in', 'ft', 'm', 'mm').
         
        :return: None.

        :raise KeyError: If the given units are not supported.
        """
        units = units.lower()
        cls.units = units_dict[units]
