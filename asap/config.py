import logging

# Initialize logger.
# Logger
with open('asap.log', 'w') as log:
    log.write('-----------------------------\n')
    log.write('ASAP LOGGING FILE INITIALIZED\n')
    log.write('-----------------------------\n')
logger = logging.getLogger('asap_logger')
logger.setLevel(logging.INFO)
_fh = logging.FileHandler('asap.log')
_fmt = logging.Formatter('%(levelname)s: %(message)s')
_fh.setFormatter(_fmt)
logger.addHandler(_fh)


class Settings(object):
    """
    Settings.

    :var float gtol: Geometric tolerance (default=1.0e-7).
    :var float ptol: Parametric tolerance (default=1.0e-9).
    :var float atol: Angular tolerance (default=1.0e-12).
    :var float stol: Shape tolerance (default=1.0e-7).
    :var float ftol: Used for curve/surface flatness criteria in
        subdivision methods (default=1.0e-3).
    :var float mtol: Mesh tolerance (default=0.005).
    :var float part_tol: Part tolerance (default=0.005).
    :var float part_angle: Part angle limit in degrees (default=30.0).
    var bool warnings: Option to print warning messages (default=*False*).
    """
    # Class variables for settings.
    gtol = 1.0e-7
    ptol = 1.0e-9
    atol = 1.0e-12
    stol = 1.0e-7
    ftol = 1.0e-3
    mtol = 0.001
    part_tol = 0.005
    part_angle = 30.0
    warnings = False

    @classmethod
    def set_gtol(cls, gtol=1.0e-7):
        """
        Set the default geometric tolerance.

        :param float gtol: Geometric tolerance.
        """
        cls.gtol = float(gtol)

    @classmethod
    def set_stol(cls, stol=1.0e-7):
        """
        Set the default shape tolerance.

        :param float stol: Shape tolerance.
        """
        cls.stol = float(stol)

    @classmethod
    def set_ptol(cls, ptol=1.0e-12):
        """
        Set the default parametric tolerance.

        :param float ptol: Parametric tolerance.
        """
        cls.ptol = float(ptol)

    @classmethod
    def set_atol(cls, atol=1.0e-12):
        """
        Set the default angular tolerance.

        :param float atol: Angular tolerance.
        """
        cls.atol = float(atol)

    @classmethod
    def set_ftol(cls, ftol=1.0e-3):
        """
        Set the default tolerance for curve/surface flatness criteria.

        :param float ftol: Flatness tolerance.
        """
        cls.ftol = float(ftol)

    @classmethod
    def set_mtol(cls, mtol=0.005):
        """
        Set the default mesh tolerance.

        :param float mtol: Mesh tolerance.
        """
        cls.mtol = float(mtol)

    @classmethod
    def set_part_tol(cls, part_tol=0.005):
        """
        Set the default frame tolerance.

        :param float part_tol: Part tolerance.
        """
        cls.part_tol = float(part_tol)

    @classmethod
    def set_part_angle(cls, part_angle=30.0):
        """
        Set the default frame angle limit.

        :param float part_angle: Part angle limit.
        """
        cls.part_angle = float(part_angle)

    @classmethod
    def activate_warnings(cls):
        """
        Activate printing warning messages.
        """
        cls.warnings = True

    @classmethod
    def deactivate_warnings(cls):
        """
        Deactivate printing warning messages.
        """
        cls.warnings = False

    @staticmethod
    def set_loggging_level(level='info'):
        """
        Set the logging level. Options include 'debug', 'info', 'warning',
        'error', or 'critical'.

        :param str level: Logging level (not case sensitive).

        :return: *True* if logging level was set, *False* if not.
        :rtype: bool

        .. note::
            Note all logging levels are currently utilized in the ASAP program.
            The default value of 'info' is recommended for general use.
        """
        log_dict = {'debug': logging.DEBUG,
                    'info': logging.INFO,
                    'warning': logging.WARNING,
                    'error': logging.ERROR,
                    'critical': logging.CRITICAL}
        if not level.lower() in log_dict:
            return False
        logger.setLevel(log_dict[level.lower()])
        return True
