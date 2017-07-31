# -*- coding: utf-8 -*-
import sys, os
from PySide.QtCore import QThread, Signal, Qt
from PySide.QtGui import QDialog, QProgressBar, QLabel, QHBoxLayout


class AutoPickThread(QThread):
    message = Signal(str)
    finished = Signal()

    def __init__(self, parent, func, infile, fnames, eventid, savepath):
        super(AutoPickThread, self).__init__()
        self.setParent(parent)
        self.func = func
        self.infile = infile
        self.fnames = fnames
        self.eventid = eventid
        self.savepath = savepath

    def run(self):
        sys.stdout = self

        picks = self.func(None, None, self.infile, self.fnames, self.eventid, self.savepath)

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
    message = Signal(str)
    
    def __init__(self, parent, func, arg=None, progressText=None, pb_widget=None, redirect_stdout=False):
        QThread.__init__(self, parent)
        self.func = func
        self.arg = arg
        self.progressText = progressText
        self.pb_widget = pb_widget
        self.redirect_stdout = redirect_stdout
        self.finished.connect(self.hideProgressbar)
        self.showProgressbar()

    def run(self):
        if self.redirect_stdout:
            sys.stdout = self        
        try:
            if self.arg:
                self.data = self.func(self.arg)
            else:
                self.data = self.func()
            self._executed = True
        except Exception as e:
            self._executed = False
            self._executedError = e
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print('Exception: {}, file: {}, line: {}'.format(exc_type, fname, exc_tb.tb_lineno))
        sys.stdout = sys.__stdout__        

    def __del__(self):
        self.wait()

    def showProgressbar(self):
        if self.progressText:
            if not self.pb_widget:
                self.pb_widget = QDialog(self.parent())
                self.pb_widget.setWindowFlags(Qt.SplashScreen)                
                self.pb_widget.setModal(True)
            hl = QHBoxLayout()
            pb = QProgressBar()
            pb.setRange(0, 0)
            hl.addWidget(pb)
            hl.addWidget(QLabel(self.progressText))
            self.pb_widget.setLayout(hl)
            self.pb_widget.show()

    def hideProgressbar(self):
        if self.pb_widget:
            self.pb_widget.hide()

    def write(self, text):
        self.message.emit(text)

    def flush(self):
        pass
    
    
class MultiThread(QThread):
    finished = Signal(str)
    message = Signal(str)    

    def __init__(self, parent, func, args, ncores=1,
                 progressText=None, pb_widget=None, redirect_stdout=False):
        QThread.__init__(self, parent)
        self.func = func
        self.args = args
        self.progressText = progressText
        self.pb_widget = pb_widget
        self.redirect_stdout = redirect_stdout
        self.finished.connect(self.hideProgressbar)
        self.showProgressbar()
        
    def run(self):
        if self.redirect_stdout:
            sys.stdout = self        
        try:
            pool = multiprocessing.Pool(self.ncores)
            self.data = pool.map_async(self.func, self.args, callback=self.emitDone)
            #self.data = pool.apply_async(self.func, self.shotlist, callback=self.emitDone) #emit each time returned
            pool.close()
            self._executed = True
            return result        
        except Exception as e:
            self._executed = False
            self._executedError = e
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            print('Exception: {}, file: {}, line: {}'.format(exc_type, fname, exc_tb.tb_lineno))
        sys.stdout = sys.__stdout__        

    def __del__(self):
        self.wait()

    def showProgressbar(self):
        if self.progressText:
            if not self.pb_widget:
                self.pb_widget = QDialog(self.parent())
                self.pb_widget.setWindowFlags(Qt.SplashScreen)                
                self.pb_widget.setModal(True)
            hl = QHBoxLayout()
            pb = QProgressBar()
            pb.setRange(0, 0)
            hl.addWidget(pb)
            hl.addWidget(QLabel(self.progressText))
            self.pb_widget.setLayout(hl)
            self.pb_widget.show()

    def hideProgressbar(self):
        if self.pb_widget:
            self.pb_widget.hide()

    def write(self, text):
        self.message.emit(text)

    def flush(self):
        pass
    
    def emitDone(self, result):
        print('emitDone!')
        self.finished.emit('Done thread!')
        self.hideProgressBar()
