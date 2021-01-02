"""
Based on https://github.com/benley/python-glog/blob/master/glog.py
(c) BSD 2-Clause 2015 Benjamin Staffin

Changelog: 2016/02  Removed gflags dependency.
"""
import logging
import os.path
import time


def format_message(record):
    try:
        record_message = '%s' % (record.msg % record.args)
    except TypeError:
        record_message = record.msg
    return record_message


class GlogFormatter(logging.Formatter):
    LEVEL_MAP = {
        logging.FATAL: 'F',
        logging.ERROR: 'E',
        logging.WARN: 'W',
        logging.INFO: 'I',
        logging.DEBUG: 'D'
    }

    def __init__(self):
        logging.Formatter.__init__(self)

    def format(self, record):
        level = GlogFormatter.LEVEL_MAP.get(record.levelno, '?')

        date = time.localtime(record.created)
        date_usec = (record.created - int(record.created)) * 1e6
        record_message = '%c%02d%02d %02d:%02d:%02d.%06d %s %s:%d] %s' % (
            level, date.tm_mon, date.tm_mday, date.tm_hour, date.tm_min,
            date.tm_sec, date_usec,
            record.process if record.process is not None else '?????',
            record.filename,
            record.lineno,
            format_message(record))
        record.getMessage = lambda: record_message
        return logging.Formatter.format(self, record)


DEBUG = logging.DEBUG
INFO = logging.INFO
WARNING = logging.WARNING
WARN = logging.WARN
ERROR = logging.ERROR
FATAL = logging.FATAL

_level_names = {
    DEBUG: 'DEBUG',
    INFO: 'INFO',
    WARN: 'WARN',
    ERROR: 'ERROR',
    FATAL: 'FATAL'
}

_level_letters = [name[0] for name in _level_names.values()]

GLOG_PREFIX_REGEX = (
                        r"""
                        (?x) ^
                        (?P<severity>[%s])
                        (?P<month>\d\d)(?P<day>\d\d)\s
                        (?P<hour>\d\d):(?P<minute>\d\d):(?P<second>\d\d)
                        \.(?P<microsecond>\d{6})\s+
                        (?P<process_id>-?\d+)\s
                        (?P<filename>[a-zA-Z<_][\w._<>-]+):(?P<line>\d+)
                        \]\s
                        """) % ''.join(_level_letters)
"""Regex you can use to parse glog line prefixes."""


def make_logger(name, level=logging.INFO):
    new_logger = logging.getLogger(name)
    # Clear any default handlers.
    if new_logger.hasHandlers():
        new_logger.handlers.clear()
    new_logger.setLevel(level)
    new_logger.propagate = False
    return new_logger


def add_file_handler(logger_obj, filename, level=logging.INFO):
    if os.path.isfile(filename):
        info('Appending to an existing log file {}'.format(filename))
    else:
        info('New log file {}'.format(filename))
    file_handler = logging.FileHandler(filename)
    file_handler.setFormatter(GlogFormatter())
    file_handler.setLevel(level)
    logger_obj.addHandler(file_handler)


def add_stream_handler(logger_obj, level=logging.INFO):
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(GlogFormatter())
    stream_handler.setLevel(level)
    logger_obj.addHandler(stream_handler)


# Clear root logger.
_root_logger = logging.getLogger()
if _root_logger.hasHandlers():
    _root_logger.handlers.clear()
_root_logger.propagate = False

logger = make_logger(name='global', level=logging.INFO)
add_stream_handler(logger, level=logging.INFO)

debug = logger.debug
info = logger.info
warning = logger.warning
warn = logger.warning
error = logger.error
exception = logger.exception
fatal = logger.fatal
log = logger.log
