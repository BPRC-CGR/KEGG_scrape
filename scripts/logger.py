import logging

# Define new log level 'RESULT'
RESULT_LEVEL_NUM = 25  # Choosing 25 as it is between INFO (20) and WARNING (30)
logging.addLevelName(RESULT_LEVEL_NUM, 'RESULT')


def result(self, message, *args, **kws):
    if self.isEnabledFor(RESULT_LEVEL_NUM):
        self._log(RESULT_LEVEL_NUM, message, args, **kws)


# Add 'result' method to Logger class
logging.Logger.result = result

# Custom formatter with colors
COLORS = {
    'WARNING': '\033[93m',
    'INFO': '\033[92m',
    'RESULT': '\033[94m',  # Light blue color for RESULT
    'DEBUG': '\033[94m',
    'CRITICAL': '\033[91m',
    'ERROR': '\033[91m',
    'ENDC': '\033[0m',
}


class CustomFormatter(logging.Formatter):
    def format(self, record):
        levelname = record.levelname
        message = logging.Formatter.format(self, record)
        return f"{COLORS.get(levelname, '')}{message}{COLORS['ENDC']}"


def custom_logger(name=__name__):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    formatter = CustomFormatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.handlers = []
    logger.addHandler(handler)
    return logger
