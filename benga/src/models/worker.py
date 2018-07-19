import threading
from concurrent.futures import ProcessPoolExecutor
from PySide import QtCore


class QSignal(QtCore.QObject):
        sig = QtCore.Signal(str)


class Worker(QtCore.QRunnable):
    def __init__(self, func, args):
        super(Worker, self).__init__()
        self.signal = QSignal()
        self.stopped = False
        self.mutex = QtCore.QMutex()
        self.func = func
        self.args = args

    def run(self, debug=False):
        if debug:
            print(threading.get_ident())
        with QtCore.QMutexLocker(self.mutex):
            self.stopped = False
        self.func(*self.args)
        # self.signal.sig.emit('OK')

    def stop(self):
        with QtCore.QMutexLocker(self.mutex):
            self.stopped = True

    def is_stopped(self):
        with QtCore.QMutexLocker(self.mutex):
            return self.stopped


class ThreadPool(QtCore.QObject):
    def __init__(self):
        self.pool = QtCore.QThreadPool.globalInstance()

    def start(self, f, args):
        self.pool.start(Worker(f, args))
        self.pool.waitForDone()
