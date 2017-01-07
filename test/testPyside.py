from PySide import QtCore, QtGui
import time


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(400, 133)
        self.progressBar = QtGui.QProgressBar(Dialog)
        self.progressBar.setGeometry(QtCore.QRect(20, 10, 361, 23))
        self.progressBar.setProperty("value", 24)
        self.progressBar.setObjectName("progressBar")
        self.pushButton = QtGui.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(20, 40, 361, 61))
        self.pushButton.setObjectName("pushButton")

        self.worker = Worker()
        self.worker.updateProgress.connect(self.setProgress)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

        self.progressBar.minimum = 1
        self.progressBar.maximum = 100

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton.setText(QtGui.QApplication.translate("Dialog", "PushButton", None, QtGui.QApplication.UnicodeUTF8))
        self.progressBar.setValue(0)
        self.pushButton.clicked.connect(self.worker.start)

    def setProgress(self, progress):
        self.progressBar.setValue(progress)


class Worker(QtCore.QThread):
    def __init__(self):
        QtCore.QThread.__init__(self)
        self.updateProgress = QtCore.Signal(int)

    def run(self):
        for i in range(1, 101):
            self.updateProgress.emit(i)
            time.sleep(0.1)

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    Dialog = QtGui.QDialog()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec_())
