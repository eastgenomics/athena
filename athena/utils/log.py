"""Functions to handle log streams"""

import logging
from sys import stdout


def get_console_handler() -> logging.StreamHandler:
    """
    Sets the stream handler to stdout with formatting

    Returns
    -------
    logging.StreamHandler
        handle for the stdout stream
    """
    console_handler = logging.StreamHandler(stdout)
    console_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(module)s] %(levelname)s: %(message)s"
        )
    )
    return console_handler


def get_logger(logger_name: str, log_level: str = None) -> logging.Logger:
    """
    Initialise the logger

    Parameters
    ----------
    logger_name : str
        name of the logger to initialise
    log_level : str
        level of logging to set

    Returns
    -------
    logging.Logger
        handle to configured logger
    """
    logger = logging.getLogger(logger_name)

    if log_level:
        logger.setLevel(log_level)

    if not logger.handlers:
        logger.addHandler(get_console_handler())

    logger.propagate = False

    return logger
