from PyQt5 import QtCore


class Worker(QtCore.QRunnable):
    def __init__(self, func, args):
        super(Worker, self).__init__()
        self.func = func
        self.args = args

    def run(self):
        self.func(*self.args)


class ThreadPool:
    def __init__(self):
        self.pool = QtCore.QThreadPool.globalInstance()

    def start(self, f, args):
        return self.pool.start(Worker(f, args))

