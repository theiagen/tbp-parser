import logging
import sys

def setup_logger(level: int = logging.INFO):
    """Configure logging for the entire application."""

    logging.basicConfig(
        encoding='utf-8',
        level=level,
        format='[%(asctime)s][%(name)s.%(funcName)s][%(levelname)s]: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.StreamHandler(sys.stderr),
            logging.FileHandler(f"tbp_parser.log", mode='a', encoding='utf-8'),
        ]
      )