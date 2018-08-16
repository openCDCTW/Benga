import logging

FMT = "(PID:%(process)d)%(asctime)-20s[%(levelname)s] %(message)s"
DATEFMT = "%Y-%m-%d %H:%M:%S"


class ConsoleLogHandler(logging.Handler):
    def __init__(self, widget):
        super(ConsoleLogHandler, self).__init__()
        self.widget = widget

    def emit(self, record):
        msg = self.format(record)
        self.widget.append(msg)

    def write(self, m):
        raise NotImplementedError()


class LoggerFactory:
    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._logger.setLevel(logging.INFO)
        self.log_format = FMT
        self.date_format = DATEFMT

    def addLogBoxHandler(self, logbox):
        logHandler = ConsoleLogHandler(logbox)
        formatter = logging.Formatter(self.log_format, self.date_format)
        logHandler.setFormatter(formatter)
        self._logger.addHandler(logHandler)

    def addFileHandler(self, logfile):
        fh = logging.FileHandler(logfile, mode="a", encoding=None, delay=False)
        formatter = logging.Formatter(self.log_format, self.date_format)
        fh.setFormatter(formatter)
        self._logger.addHandler(fh)

    def addConsoleHandler(self):
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter(self.log_format, self.date_format)
        ch.setFormatter(formatter)
        self._logger.addHandler(ch)

    def redircet_stdout(self):
        pass

    def create(self):
        return self._logger

