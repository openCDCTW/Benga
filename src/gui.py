import sys
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from utils import *
import algorithms


def set_enable(widgets, enabled):
    for widget in widgets.values():
        widget.setEnabled(enabled)


class Window(QWidget):
    def __init__(self):
        super(Window, self).__init__()
        self.createMainLayout()

    def createMainLayout(self):
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
        self.inputPanel1 = {}
        self.inputPanel1["label"] = QLabel("Directory:")
        self.inputPanel1["label"].setFixedSize(80, 30)
        self.inputPanel1["line"] = QLineEdit()
        self.inputPanel1["browseButton"] = QPushButton("&Browse...")
        self.inputPanel1["browseButton"].setFixedSize(80, 30)
        self.inputPanel1["browseButton"].clicked.connect(lambda: self.browse(self.inputPanel1))

        layout = QHBoxLayout()
        layout.addWidget(self.inputPanel1["label"])
        layout.addWidget(self.inputPanel1["line"])
        layout.addWidget(self.inputPanel1["browseButton"])
        return layout

    def createInput2Layout(self):
        self.inputPanel2 = {}
        self.inputPanel2["label"] = QLabel("Directory:")
        self.inputPanel2["label"].setFixedSize(80, 30)
        self.inputPanel2["line"] = QLineEdit()
        self.inputPanel2["browseButton"] = QPushButton("&Browse...")
        self.inputPanel2["browseButton"].setFixedSize(80, 30)
        self.inputPanel2["browseButton"].clicked.connect(lambda: self.browse(self.inputPanel2))

        layout = QHBoxLayout()
        layout.addWidget(self.inputPanel2["label"])
        layout.addWidget(self.inputPanel2["line"])
        layout.addWidget(self.inputPanel2["browseButton"])
        return layout

    def createSubmit1Layout(self):
        self.submitButton1 = QPushButton("&Submit")
        self.submitButton1.setFixedSize(80, 30)
        self.submitButton1.clicked.connect(lambda: self.submit(self.inputPanel1))

        layout = QHBoxLayout()
        layout.addWidget(self.submitButton1)
        layout.setAlignment(Qt.AlignRight)
        return layout

    def createSubmit2Layout(self):
        self.submitButton2 = QPushButton("&Submit")
        self.submitButton2.setFixedSize(80, 30)
        self.submitButton2.clicked.connect(lambda: self.submit(self.inputPanel2))

        layout = QHBoxLayout()
        layout.addWidget(self.submitButton2)
        layout.setAlignment(Qt.AlignRight)
        return layout

    def browse(self, panel):
        dialog = QFileDialog()
        dialog.setFileMode(QFileDialog.Directory)
        dialog.setOption(QFileDialog.ShowDirsOnly, True)
        directory = dialog.getExistingDirectory(dialog, "Find Directory", QDir.currentPath())

        if directory:
            panel["line"].setText(directory)

    def submit(self, panel):
        source_dir = str(panel["line"].text())
        parent_dir = os.path.dirname(source_dir)
        working_dir = os.path.join(parent_dir, "temp")
        clear_folder(working_dir)
        set_enable(panel, False)
        algorithms.main(working_dir, source_dir)
        set_enable(panel, True)


if __name__ == "__main__":
    app = QApplication(sys.argv)

    w = Window()
    w.resize(500, 150)
    w.setWindowTitle("MSGA")
    w.show()

    sys.exit(app.exec_())
