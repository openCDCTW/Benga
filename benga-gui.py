from PySide import QtGui, QtCore
from src.algorithms import databases, phylogeny, profiling
from src.models import worker
from src.models.jobs import *
from src.models.logs import *
from src.utils import files
from src.view import Ui_mainTab


class Window(Ui_mainTab):
    def __init__(self, mainTab):
        super(Window, self).__init__()
        self.setupUi(mainTab)

        current_dir = os.path.dirname(__file__)
        self.ROOT_DIR = os.path.abspath(os.path.join(current_dir, os.pardir))
        self.data_dir = os.path.join(self.ROOT_DIR, "data")
        files.create_if_not_exist(self.data_dir)

        self.pool = worker.ThreadPool()

        # setting behaviers
        self.fileSelector.clicked.connect(lambda: self.select_dir(self.fileSelectText))
        self.fileSubmitter.clicked.connect(lambda: self.pool.start(self.make_database, ()))

        self.profileSelector.clicked.connect(lambda: self.select_profiles(self.profileSelectText, self.joblist_2))
        self.runButton.clicked.connect(lambda: self.pool.start(self.make_profile, ()))

        self.plottingSelector.clicked.connect(lambda: self.select_dir(self.plottingSelectText))
        self.plottingSelectText.textChanged.connect(lambda: self.plotDendrogram())

    def select_dir(self, text_widget):
        dialog = QtGui.QFileDialog()
        dialog.setFileMode(QtGui.QFileDialog.Directory)
        dialog.setOption(QtGui.QFileDialog.ShowDirsOnly, True)
        directory = dialog.getExistingDirectory(dialog, "Select Directory", QtCore.QDir.currentPath())
        if directory:
            text_widget.setText(directory)

    def show_dbs(self, db_box):
        database_dir = files.joinpath(self.data_dir, JobType.PGDB.to_str())
        dbs = [x for x in os.listdir(database_dir) if os.path.isdir(files.joinpath(database_dir, x))]
        for i in dbs:
            item = QtGui.QListWidgetItem(i)
            db_box.addItem(item)
        db_box.show()

    def select_profiles(self, text_widget, job_box):
        self.select_dir(text_widget)
        self.show_dbs(job_box)

    def make_database(self):
        """take files from fileSelectText for making database"""
        switch_widgets = [self.fileSubmitter, self.fileSelectText, self.fileSelector]
        self.disable(switch_widgets)

        # setup paths
        source_dir = str(self.fileSelectText.toPlainText())
        database_dir = files.joinpath(self.data_dir, JobType.PGDB.to_str())
        files.create_if_not_exist(database_dir)

        # create new job
        self.jobmgr = JobManager(self.data_dir)
        jobid, job_dir = self.jobmgr.start_job(JobType.PGDB)
        self.jobmgr.close()

        # setup logger
        factory = LoggerFactory()
        factory.addLogBoxHandler(self.logbox_1)  # TODO: bug -- crash after Worker is done.
        factory.addFileHandler(files.joinpath(database_dir, "log_" + jobid + ".txt"))
        logger = factory.create()

        # process algorithms
        databases.annotate_configs(source_dir, job_dir, logger=logger)
        databases.make_database(job_dir, logger=logger)

        self.enable(switch_widgets)

    def make_profile(self):
        """take file from profileSelectText for profiling"""
        switch_widgets = [self.profileSelector, self.profileSelectText, self.joblist_2, self.runButton]
        self.disable(switch_widgets)

        # setup paths
        option = str(self.joblist_2.currentItem().text())
        database_dir = files.joinpath(self.data_dir, JobType.PGDB.to_str(), option, "DB")
        profile_dir = str(self.profileSelectText.toPlainText())
        query_dir = files.joinpath(self.data_dir, JobType.WGMLST.to_str())
        files.create_if_not_exist(query_dir)

        # create new job
        self.jobmgr = JobManager(self.data_dir)
        jobid, job_dir = self.jobmgr.start_job(JobType.WGMLST)
        self.jobmgr.close()

        # setup logger
        factory = LoggerFactory()
        factory.addLogBoxHandler(self.logbox_2)  # TODO: bug -- crash after Worker is done.
        factory.addFileHandler(files.joinpath(query_dir, "log_" + jobid + ".txt"))
        logger = factory.create()

        # process algorithms
        profiling.profiling(job_dir, profile_dir, database_dir, logger)

        self.enable(switch_widgets)

    def plotDendrogram(self):
        """take file from plottingSelectText, plot dendrogram to treeViewer"""
        switch_widgets = [self.plottingSelector, self.plottingSelectText]
        self.disable(switch_widgets)

        profile_dir = str(self.plottingSelectText.toPlainText())
        self.tree, self.newick = phylogeny.make_tree(profile_dir, "wgMLST_pan")
        self.tree.show()  # TODO: bug gui fixed

        self.enable(switch_widgets)

    @classmethod
    def enable(cls, widgets):
        for w in widgets:
            w.setEnabled(True)

    @classmethod
    def disable(cls, widgets):
        for w in widgets:
            w.setEnabled(False)


if __name__ == "__main__":
    import sys
    import logging
    logging.basicConfig()
    app = QtGui.QApplication(sys.argv)
    mainTab = QtGui.QTabWidget()
    ui = Window(mainTab)
    mainTab.show()
    sys.exit(app.exec_())
