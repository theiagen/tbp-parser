import logging
import sys
from pathlib import Path

def setup_logger(
    output_prefix: Path,
    level: int,
) -> None:
    """Configure logging for the entire application."""

    logging.basicConfig(
        encoding='utf-8',
        level=level,
        format='[%(asctime)s][%(name)s.%(funcName)s][%(levelname)s]: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler(f"{output_prefix}.log", mode='w', encoding='utf-8'),
        ]
      )