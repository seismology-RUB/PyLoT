# -*- coding: utf-8 -*-
import sys
from PySide.QtCore import QThread, Signal, Qt
from PySide.QtGui import QDialog, QProgressBar, QLabel, QVBoxLayout


class AutoPickThread(QThread):
    message = Signal(str)
    finished = Signal()

    def __init__(self, parent, func, infile, fnames, savepath):
        super(AutoPickThread, self).__init__()
        self.setParent(parent)
        self.func = func
        self.infile = infile
        self.fnames = fnames
        self.savepath = savepath

    def run(self):
        sys.stdout = self

        picks = self.func(self.infile, self.fnames, self.savepath)

        print("Autopicking finished!\n")

        try:
            for station in picks:
                self.parent().addPicks(station, picks[station], type='auto')
        except AttributeError:
            print(picks)
        sys.stdout = sys.__stdout__
        self.finished.emit()

    def write(self, text):
        self.message.emit(text)

    def flush(self):
        pass


class Thread(QThread):
    def __init__(self, parent, func, progressText = None):
        QThread.__init__(self, parent)
        self.func = func
        self.progressText = progressText
        self.pbdlg = None
        self.finished.connect(self.hideProgressbar)
        self.showProgressbar()

    def run(self):
        self.func()

    def __del__(self):
        self.wait()

    def showProgressbar(self):
        if self.progressText:
            self.pbdlg = QDialog(self.parent())
            self.pbdlg.setModal(True)
            vl = QVBoxLayout()
            pb = QProgressBar()
            pb.setRange(0, 0)
            vl.addWidget(pb)
            vl.addWidget(QLabel(self.progressText))
            self.pbdlg.setLayout(vl)
            self.pbdlg.show()
            self.pbdlg.setWindowFlags(Qt.FramelessWindowHint)
            self.pbdlg.show()

    def hideProgressbar(self):
        if self.pbdlg:
            self.pbdlg.hide()

            
