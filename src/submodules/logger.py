
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(filename)s %(funcName)s:%(lineno)d - %(levelname)s: %(message)s'
)
LOGGER = logging.getLogger()


def quiet_log_info(quiet, *args):
    if not quiet:
        LOGGER.info(*args, stacklevel=2)


def quiet_log_warning(quiet, *args):
    if not quiet:
        LOGGER.warning(*args, stacklevel=2)


def quiet_log_error(quiet, *args):
    if not quiet:
        LOGGER.error(*args, stacklevel=2)
