import logging

# Initialize logger.
with open('afem.log', 'w') as log:
    log.write('-----------------------------\n')
    log.write('AFEM LOGGING FILE INITIALIZED\n')
    log.write('-----------------------------\n')
logger = logging.getLogger('afem_logger')
logger.setLevel(logging.INFO)
_fh = logging.FileHandler('afem.log')
_fmt = logging.Formatter('%(levelname)s: %(message)s')
_fh.setFormatter(_fmt)
logger.addHandler(_fh)
logger.propagate = False

# Dictionary for units
units_dict = {'i': 'IN',
              'in': 'IN',
              'inch': 'IN',
              'inches': 'IN',
              'f': 'FT',
              'ft': 'FT',
              'foot': 'FT',
              'feet': 'FT',
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
    units = 'IN'

    @staticmethod
    def log_to_console(option=False):
        """
        Option to print logging output to the console.

        :param bool option: *True* will print to the console, *False* will
            not.

        :return: None.
        """
        if option:
            logger.propagate = True
        else:
            logger.propagate = False

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
