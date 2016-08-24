import logging
import sys

# Get root logger (all other loggers will be derived from this logger's
# properties)
logger = logging.getLogger()
logger.handlers[0].stream = sys.stdout
