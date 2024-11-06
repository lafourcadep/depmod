from __future__ import annotations

import logging
import sys

def get_logger(level: int = 10) -> logging.Logger:
    logger = logging.getLogger(__name__)
    


    logger.setLevel(level)
    logger.propagate = False

    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setLevel(level)
    consoleHandler.setFormatter(
        logging.Formatter(
            "[%(levelname)s] %(message)s"
        )
    )

    logger.addHandler(consoleHandler)

    return logger

logger = get_logger(10)
