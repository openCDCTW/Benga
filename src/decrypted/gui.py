import os
from PyQt5 import QtWidgets, QtCore
from design import Ui_mainTab
from algorithms import *
from utils import *
from jobs import *
from logs import *
from processes import *


class Window(Ui_mainTab):
    def __init__(self, mainTab):
        super(Window, self).__init__()
        self.setupUi(mainTab)

        current_dir = os.path.dirname(__file__)
        self.ROOT_DIR = os.path.abspath(os.path.join(current_dir, os.pardir))
        self.data_dir = os.path.join(self.ROOT_DIR, "data")
        create_if_not_exist(self.data_dir)

        # initialize managers
        self.pool = ThreadPool()

        # setting behaviers
        self.fileSelector.clicked.connect(lambda: self.select_dir(self.fileSelectText))
        self.fileSubmitter.clicked.connect(lambda: self.pool.start(self.submitForDB, ()))

        self.profileSelector.clicked.connect(lambda: self.select_profiles(self.profileSelectText, self.joblist_2))
        self.runButton.clicked.connect(lambda: self.pool.start(self.submitForProfile, ()))

        self.plottingSelector.clicked.connect(lambda: self.select_dir(self.plottingSelectText))
        self.plottingSelectText.textChanged.connect(lambda: self.plotDendrogram())
        self.toPdfButton.clicked.connect(lambda: self.save2pdf())
        self.toNewickButton.clicked.connect(lambda: self.save2newick())

    def select_dir(self, text_widget):
        dialog = QtWidgets.QFileDialog()
        dialog.setFileMode(QtWidgets.QFileDialog.Directory)
        dialog.setOption(QtWidgets.QFileDialog.ShowDirsOnly, True)
        directory = dialog.getExistingDirectory(dialog, "Select Directory", QtCore.QDir.currentPath())
        if directory:
            text_widget.setText(directory)

    def show_dbs(self, db_box):
        database_dir = joinpath(self.data_dir, JobType.PGDB.to_str())
        dbs = [x for x in os.listdir(database_dir) if os.path.isdir(joinpath(database_dir, x))]
        for i in dbs:
            item = QtWidgets.QListWidgetItem(i)
            db_box.addItem(item)
        db_box.show()

    def select_profiles(self, text_widget, job_box):
        self.select_dir(text_widget)
        self.show_dbs(job_box)

    def submitForDB(self):
        """take files from fileSelectText for making database"""
        switch_widgets = [self.fileSubmitter, self.fileSelectText, self.fileSelector]
        self.set_enable(switch_widgets, False)

        # setup paths
        source_dir = str(self.fileSelectText.toPlainText())
        database_dir = joinpath(self.data_dir, JobType.PGDB.to_str())
        create_if_not_exist(database_dir)

        # create new job
        self.jobmgr = JobManager()
        # jobid, job_dir = self.jobmgr.start_job(JobType.PGDB)
        jobid = "d7f01d24-c0d4-11e6-8d68-dc85de763f2b" #self.jobmgr.jobs[0].jobid
        job_dir = joinpath(database_dir, jobid)
        self.jobmgr.close()

        # setup logger
        factory = LoggerFactory()
        factory.addLogBoxHandler(self.logbox_1)  # TODO: bug -- crash after Worker is done.
        factory.addFileHandler(joinpath(database_dir, "log_" + jobid + ".txt"))
        logger = factory.create()

        # process algorithms
        annotate_configs(source_dir, job_dir, logger=logger)
        make_profiles(job_dir, logger=logger)

        self.set_enable(switch_widgets, True)

    def submitForProfile(self):
        """take file from profileSelectText for profiling"""
        switch_widgets = [self.profileSelector, self.profileSelectText, self.joblist_2, self.runButton]
        self.set_enable(switch_widgets, False)

        # setup paths
        option = str(self.joblist_2.currentItem().text())
        database_dir = joinpath(self.data_dir, JobType.PGDB.to_str(), option, "DB")
        profile_dir = str(self.profileSelectText.toPlainText())
        query_dir = joinpath(self.data_dir, JobType.WGMLST.to_str())
        create_if_not_exist(query_dir)

        # create new job
        self.jobmgr = JobManager()
        jobid, job_dir = self.jobmgr.start_job(JobType.WGMLST)
        # jobid = self.jobmgr.jobs[1].jobid
        # job_dir = joinpath(query_dir, jobid)
        self.jobmgr.close()

        # setup logger
        factory = LoggerFactory()
        factory.addLogBoxHandler(self.logbox_2)  # TODO: bug -- crash after Worker is done.
        factory.addFileHandler(joinpath(query_dir, "log_" + jobid + ".txt"))
        logger = factory.create()

        # process algorithms
        wg_profiling(job_dir, profile_dir, database_dir, logger)

        self.set_enable(switch_widgets, True)

    def plotDendrogram(self):
        """take file from plottingSelectText, plot dendrogram to treeViewer"""
        switch_widgets = [self.plottingSelector, self.plottingSelectText, self.toPdfButton, self.toNewickButton]
        self.set_enable(switch_widgets, False)

        profile_dir = str(self.plottingSelectText.toPlainText())

        self.tree, self.newick = phylogenetic_tree(profile_dir, "wgMLST_core")
        self.tree.show()
        # self.treeViewer

        self.set_enable(switch_widgets, True)

    def save2pdf(self):
        pass

    def save2newick(self):
        pass

    @classmethod
    def set_enable(cls, widgets, option):
        for w in widgets:
            w.setEnabled(option)


if __name__ == "__main__":
    import sys
    logging.basicConfig()
    app = QtWidgets.QApplication(sys.argv)
    mainTab = QtWidgets.QTabWidget()
    ui = Window(mainTab)
    mainTab.show()
    sys.exit(app.exec_())
