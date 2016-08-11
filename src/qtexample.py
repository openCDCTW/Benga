import sys
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from utils import *
import contigAnnotation


class Window(QWidget):
    def __init__(self):
        super(Window, self).__init__()
        self.createLayout()

    def createLayout(self):

        layout = QVBoxLayout()
        layout.addLayout(self.createTitle1Layout())
        layout.addLayout(self.createInput1Layout())
        layout.addLayout(self.createSubmit1Layout())
        layout.addLayout(self.createTitle2Layout())
        layout.addLayout(self.createInput2Layout())
        layout.addLayout(self.createSubmit2Layout())
        self.setLayout(layout)

    def createTitle1Layout(self):
        self.title1 = QLabel("Make database")
        self.title1.setFixedSize(120, 30)
        layout = QHBoxLayout()
        layout.addWidget(self.title1)
        layout.setAlignment(Qt.AlignLeft)
        return layout

    def createTitle2Layout(self):
        self.title2 = QLabel("Build tree")
        self.title2.setFixedSize(120, 30)
        layout = QHBoxLayout()
        layout.addWidget(self.title2)
        layout.setAlignment(Qt.AlignLeft)
        return layout

    def createInput1Layout(self):
        self.label1 = QLabel("Directory:")
        self.label1.setFixedSize(80, 30)
        self.line1 = QLineEdit()
        self.browseButton1 = QPushButton("&Browse...")
        self.browseButton1.setFixedSize(80, 30)
        self.browseButton1.clicked.connect(lambda: self.browse(self.line1))

        layout = QHBoxLayout()
        layout.addWidget(self.label1)
        layout.addWidget(self.line1)
        layout.addWidget(self.browseButton1)
        return layout

    def createInput2Layout(self):
        self.label2 = QLabel("Directory:")
        self.label2.setFixedSize(80, 30)
        self.line2 = QLineEdit()
        self.browseButton2 = QPushButton("&Browse...")
        self.browseButton2.setFixedSize(80, 30)
        self.browseButton2.clicked.connect(lambda: self.browse(self.line2))

        layout = QHBoxLayout()
        layout.addWidget(self.label2)
        layout.addWidget(self.line2)
        layout.addWidget(self.browseButton2)
        return layout

    def createSubmit1Layout(self):
        self.submitButton1 = QPushButton("&Submit")
        self.submitButton1.setFixedSize(80, 30)
        self.submitButton1.clicked.connect(lambda: self.submit(self.line1))

        layout = QHBoxLayout()
        layout.addWidget(self.submitButton1)
        layout.setAlignment(Qt.AlignRight)
        return layout

    def createSubmit2Layout(self):
        self.submitButton2 = QPushButton("&Submit")
        self.submitButton2.setFixedSize(80, 30)
        self.submitButton2.clicked.connect(lambda: self.submit(self.line2))

        layout = QHBoxLayout()
        layout.addWidget(self.submitButton2)
        layout.setAlignment(Qt.AlignRight)
        return layout

    def browse(self, line):
        dialog = QFileDialog()
        dialog.setFileMode(QFileDialog.Directory)
        dialog.setOption(QFileDialog.ShowDirsOnly, True)
        directory = dialog.getExistingDirectory(dialog, "Find Directory", QDir.currentPath())

        if directory:
            line.setText(directory)

    def submit(self, line):
        source_dir = str(line.text())
        parent_dir = os.path.dirname(source_dir)
        working_dir = os.path.join(parent_dir, "temp")
        clear_folder(working_dir)
        self.line1.setEnabled(False)
        self.submitButton1.setEnabled(False)
        contigAnnotation.main(working_dir, source_dir)

    @classmethod
    def enable(cls, widgets, enabled):
        for w in widgets.values():
            w.setEnabled(enabled)


if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = Window()
    w.resize(500, 150)
    w.setWindowTitle("MSGA")
    w.show()

    sys.exit(app.exec_())
