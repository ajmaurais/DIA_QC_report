import sys
from datetime import datetime
import logging
from typing import Any
from logging import DEBUG, INFO, WARNING, ERROR, CRITICAL

_FORMATING_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR']


class _SplitStreamHandler(logging.Handler):
    '''
    Send < ERROR to stdout and >= ERROR to stderr with a custom format that
    references the parent Logger's live configuration (so toggling options at
    runtime is reflected immediately).
    '''
    def __init__(self, parent: "Logger") -> None:
        super().__init__(level=parent.level)
        self._parent = parent


    def emit(self, record: logging.LogRecord) -> None:
        if self._parent.quiet:
            return

        parts: list[str] = []
        if self._parent.show_date:
            parts.append(datetime.fromtimestamp(record.created).strftime(self._parent.dateformat))
            parts.append('-')

        if self._parent.debug_mode:
            parts.append(f"{record.filename} {record.funcName}:{record.lineno} ")

        if self._parent.show_level:
            if len(parts) > 0:
                parts.append(f'[{record.levelname}]')
            else:
                parts.append(f'[{record.levelname}]'.ljust(max(len(level) for level in _FORMATING_LEVELS)))

        parts.append(record.getMessage())

        stream = sys.stderr if record.levelno >= logging.ERROR else sys.stdout
        stream.write(" ".join(parts) + "\n")
        stream.flush()


class Logger(logging.Logger):
    '''
    Custom logger class built on top of the standard logging.Logger.
    '''
    def __init__(
        self,
        name: str = "project",
        *,
        dateformat: str = "%Y-%m-%d %H:%M:%S",
        min_level: int | str = logging.INFO,
        debug: bool = False,
        quiet: bool = False,
        show_level: bool = True,
        show_date: bool = False,
    ) -> None:
        super().__init__(name, level=min_level)

        self.dateformat = dateformat
        self.debug_mode = debug
        self.quiet = quiet
        self.show_level = show_level
        self.show_date = show_date
        self._install_handler()

        self.propagate = False


    def _install_handler(self) -> None:
        '''Remove any existing custom handlers and install a fresh one.'''
        for h in list(self.handlers):
            self.removeHandler(h)
        self.addHandler(_SplitStreamHandler(self))

    def set_debug(self, flag: bool = True) -> None:
        self.debug_mode = flag

    def set_quiet(self, flag: bool = True) -> None:
        self.quiet = flag

    def set_min_level(self, level: int | str) -> None:
        self.setLevel(level)

    info = logging.Logger.info
    debug = logging.Logger.debug
    warning = logging.Logger.warning
    warn = logging.Logger.warning
    error = logging.Logger.error


logging.setLoggerClass(Logger)

LOGGER = logging.getLogger('DIA_QC_report')


def quiet_log_info(quiet, *args):
    if not quiet:
        LOGGER.info(*args, stacklevel=2)


def quiet_log_warning(quiet, *args):
    if not quiet:
        LOGGER.warning(*args, stacklevel=2)


def quiet_log_error(quiet, *args):
    if not quiet:
        LOGGER.error(*args, stacklevel=2)