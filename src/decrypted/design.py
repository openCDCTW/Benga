# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'gui/design.ui'
#
# Created by: PyQt5 UI code generator 5.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_mainTab(object):
    def setupUi(self, mainTab):
        mainTab.setObjectName("mainTab")
        mainTab.resize(534, 424)
        self.panGenomeTab = QtWidgets.QWidget()
        self.panGenomeTab.setObjectName("panGenomeTab")
        self.fileSelectText = QtWidgets.QTextEdit(self.panGenomeTab)
        self.fileSelectText.setGeometry(QtCore.QRect(30, 30, 361, 31))
        self.fileSelectText.setObjectName("fileSelectText")
        self.fileSelector = QtWidgets.QPushButton(self.panGenomeTab)
        self.fileSelector.setGeometry(QtCore.QRect(409, 30, 101, 29))
        self.fileSelector.setObjectName("fileSelector")
        self.logbox_1 = QtWidgets.QTextBrowser(self.panGenomeTab)
        self.logbox_1.setGeometry(QtCore.QRect(30, 80, 361, 281))
        self.logbox_1.setObjectName("messagebox_1")
        self.fileSubmitter = QtWidgets.QPushButton(self.panGenomeTab)
        self.fileSubmitter.setGeometry(QtCore.QRect(410, 80, 101, 29))
        self.fileSubmitter.setObjectName("fileSubmitter")
        self.jobFetcher = QtWidgets.QPushButton(self.panGenomeTab)
        self.jobFetcher.setGeometry(QtCore.QRect(409, 330, 101, 29))
        self.jobFetcher.setObjectName("jobFetcher")
        self.joblist_1 = QtWidgets.QListWidget(self.panGenomeTab)
        self.joblist_1.setGeometry(QtCore.QRect(410, 120, 101, 201))
        self.joblist_1.setObjectName("joblist_1")
        mainTab.addTab(self.panGenomeTab, "")

        self.profilerTab = QtWidgets.QWidget()
        self.profilerTab.setObjectName("profilerTab")
        self.profileSelectText = QtWidgets.QTextEdit(self.profilerTab)
        self.profileSelectText.setGeometry(QtCore.QRect(30, 30, 361, 31))
        self.profileSelectText.setObjectName("profileSelectText")
        self.profileSelector = QtWidgets.QPushButton(self.profilerTab)
        self.profileSelector.setGeometry(QtCore.QRect(409, 30, 101, 29))
        self.profileSelector.setObjectName("profileSelector")
        self.joblist_2 = QtWidgets.QListWidget(self.profilerTab)
        self.joblist_2.setGeometry(QtCore.QRect(220, 80, 171, 281))
        self.joblist_2.setObjectName("joblist_2")
        self.logbox_2 = QtWidgets.QTextBrowser(self.profilerTab)
        self.logbox_2.setGeometry(QtCore.QRect(30, 80, 171, 281))
        self.logbox_2.setObjectName("messagebox_2")
        self.runButton = QtWidgets.QPushButton(self.profilerTab)
        self.runButton.setGeometry(QtCore.QRect(410, 330, 101, 29))
        self.runButton.setObjectName("runButton")
        mainTab.addTab(self.profilerTab, "")

        self.dendrogramTab = QtWidgets.QWidget()
        self.dendrogramTab.setObjectName("dendrogramTab")
        self.treeViewer = QtWidgets.QGraphicsView(self.dendrogramTab)
        self.treeViewer.setGeometry(QtCore.QRect(30, 80, 361, 291))
        self.treeViewer.setObjectName("treeViewer")
        self.plottingSelector = QtWidgets.QPushButton(self.dendrogramTab)
        self.plottingSelector.setGeometry(QtCore.QRect(410, 30, 101, 29))
        self.plottingSelector.setObjectName("plottingSelector")
        self.plottingSelectText = QtWidgets.QTextEdit(self.dendrogramTab)
        self.plottingSelectText.setGeometry(QtCore.QRect(30, 30, 361, 31))
        self.plottingSelectText.setObjectName("plottingSelectText")
        self.toPdfButton = QtWidgets.QPushButton(self.dendrogramTab)
        self.toPdfButton.setGeometry(QtCore.QRect(410, 290, 101, 29))
        self.toPdfButton.setObjectName("toPdfButton")
        self.toNewickButton = QtWidgets.QPushButton(self.dendrogramTab)
        self.toNewickButton.setGeometry(QtCore.QRect(410, 340, 101, 29))
        self.toNewickButton.setObjectName("toNewickButton")
        mainTab.addTab(self.dendrogramTab, "")

        self.retranslateUi(mainTab)
        mainTab.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(mainTab)
        mainTab.setTabOrder(self.fileSelectText, self.fileSelector)
        mainTab.setTabOrder(self.fileSelector, self.logbox_1)
        mainTab.setTabOrder(self.logbox_1, self.fileSubmitter)
        mainTab.setTabOrder(self.fileSubmitter, self.joblist_1)
        mainTab.setTabOrder(self.joblist_1, self.jobFetcher)
        mainTab.setTabOrder(self.jobFetcher, self.profileSelectText)
        mainTab.setTabOrder(self.profileSelectText, self.profileSelector)
        mainTab.setTabOrder(self.profileSelector, self.logbox_2)
        mainTab.setTabOrder(self.logbox_2, self.joblist_2)
        mainTab.setTabOrder(self.joblist_2, self.runButton)
        mainTab.setTabOrder(self.runButton, self.plottingSelectText)
        mainTab.setTabOrder(self.plottingSelectText, self.plottingSelector)
        mainTab.setTabOrder(self.plottingSelector, self.treeViewer)
        mainTab.setTabOrder(self.treeViewer, self.toPdfButton)
        mainTab.setTabOrder(self.toPdfButton, self.toNewickButton)

    def retranslateUi(self, mainTab):
        _translate = QtCore.QCoreApplication.translate
        mainTab.setWindowTitle(_translate("mainTab", "wgMLST GUI"))
        self.fileSelectText.setPlaceholderText(_translate("mainTab", "select files for buidling database..."))
        self.fileSelector.setText(_translate("mainTab", "Select Files"))
        self.fileSubmitter.setText(_translate("mainTab", "Submit"))
        self.jobFetcher.setText(_translate("mainTab", "Fetch Job"))
        mainTab.setTabText(mainTab.indexOf(self.panGenomeTab), _translate("mainTab", "Pan-genome DB"))
        self.profileSelectText.setPlaceholderText(_translate("mainTab", "select profile..."))
        self.profileSelector.setText(_translate("mainTab", "Select Profile"))
        self.runButton.setText(_translate("mainTab", "Run"))
        mainTab.setTabText(mainTab.indexOf(self.profilerTab), _translate("mainTab", "wgMLST Profiler"))
        self.plottingSelector.setText(_translate("mainTab", "Select Profile"))
        self.plottingSelectText.setPlaceholderText(_translate("mainTab", "select profile for plotting..."))
        self.toPdfButton.setText(_translate("mainTab", "Save as pdf"))
        self.toNewickButton.setText(_translate("mainTab", "Save as newick"))
        mainTab.setTabText(mainTab.indexOf(self.dendrogramTab), _translate("mainTab", "Dendrogram"))
