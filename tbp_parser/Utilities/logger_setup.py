import logging
import sys
from pathlib import Path

def setup_logger(
    log_file: Path | None = None,
    level: int = logging.INFO,
) -> None:
    """Configure logging for the entire application. Called once, from main().

    Args:
        log_file (Path | None): the full path of the log file to write; if None, only stderr is used
        level (int): the logging level to emit at
    """
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file is not None:
        handlers.append(logging.FileHandler(f"{log_file}.log", mode='w', encoding='utf-8'))

    logging.basicConfig(
        level=level,
        format='[%(asctime)s][%(name)s.%(funcName)s][%(levelname)s]: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=handlers,
      )