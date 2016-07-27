#!/usr/bin/env python
# -*- coding: utf-8 -*-
from PySide import QtCore, QtGui
from pylot.core.active import surveyUtils, activeSeismoPick, seismicArrayPreparation
from generate_survey_layout import *
from generate_survey_layout_minimal import *
from generate_seisarray_layout import *
from picking_parameters_layout import *

import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

def openFile(name = 'Open'):
    dialog = QtGui.QFileDialog()
    dialog.setWindowTitle(name)                #not working yet
    filename = dialog.getOpenFileName()
    if len(filename[0]) > 0:
        return filename[0]

def saveFile(name = 'Save'):
    dialog = QtGui.QFileDialog()
    dialog.setWindowTitle(name)
    filename = dialog.getSaveFileName()
    if len(filename[0]) > 0:
        return filename[0]

def browseDir(name = 'Open Directory'):
    dialog = QtGui.QFileDialog()
    dialog.setWindowTitle(name)
    directory = dialog.getExistingDirectory()
    if len(directory) > 0:
        return directory

def getMaxCPU():
    import multiprocessing
    return multiprocessing.cpu_count()

class Gen_SeisArray(object):
    def __init__(self, mainwindow):
        self.mainwindow = mainwindow
        self.seisarray = None
        self.srcfile = None
        self.recfile = None
        self.ptsfile = None
        self.init_dialog()
        self.start_dialog()

    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_generate_seisarray()
        ui.setupUi(qdialog)
        self.ui = ui
        self.qdialog = qdialog
        self.connectButtons()

    def start_dialog(self):
        self.init_last_selection()
        if self.qdialog.exec_():
            self.refresh_filenames()
            if self.ui.radioButton_interpolatable.isChecked():
                self.seisarray = seismicArrayPreparation.SeisArray(self.recfile, True)
            elif self.ui.radioButton_normal.isChecked():
                self.seisarray = seismicArrayPreparation.SeisArray(self.recfile, False)
            if len(self.srcfile) > 0:
                self.seisarray.addSourceLocations(self.srcfile)
            if len(self.ptsfile) > 0:
                self.seisarray.addMeasuredTopographyPoints(self.ptsfile)
            self.executed = True
        else:
            self.refresh_filenames()
            self.executed = False

    def refresh_filenames(self):
        self.srcfile = self.ui.lineEdit_src.text()
        self.recfile = self.ui.lineEdit_rec.text()
        self.ptsfile = self.ui.lineEdit_pts.text()

    def init_last_selection(self):
        self.ui.lineEdit_src.setText(self.srcfile)
        self.ui.lineEdit_rec.setText(self.recfile)
        self.ui.lineEdit_pts.setText(self.ptsfile)

    def get_seisarray(self):
        if self.seisarray is not None:
            return self.seisarray

    def connectButtons(self):
        QtCore.QObject.connect(self.ui.pushButton_rec, QtCore.SIGNAL("clicked()"), self.chooseMeasuredRec)
        QtCore.QObject.connect(self.ui.pushButton_src, QtCore.SIGNAL("clicked()"), self.chooseMeasuredSrc)
        QtCore.QObject.connect(self.ui.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseMeasuredPts)

    def chooseMeasuredSrc(self):
        self.ui.lineEdit_src.setText(openFile('Open measured sources file.'))

    def chooseMeasuredRec(self):
        self.ui.lineEdit_rec.setText(openFile('Open measured receivers file.'))

    def chooseMeasuredPts(self):
        self.ui.lineEdit_pts.setText(openFile('Open measured points file.'))



class Gen_Survey_from_SA(object):
    def __init__(self, mainwindow, seisarray):
        self.mainwindow = mainwindow
        self.seisarray = seisarray
        self.survey = None
        self.obsdir = None
        self.fstart = 'shot'
        self.fend = '.dat'
        self.init_dialog()
        self.start_dialog()

    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_generate_survey_minimal()
        ui.setupUi(qdialog)
        self.ui = ui
        self.qdialog = qdialog
        self.connectButtons()

    def start_dialog(self):
        self.init_last_selection()
        if self.qdialog.exec_():
            self.refresh_filenames()
            self.survey = activeSeismoPick.Survey(self.obsdir, seisArray = self.seisarray,
                                                  useDefaultParas = True, fstart = self.fstart,
                                                  fend = self.fend)
            self.executed = True
        else:
            self.refresh_filenames()
            self.executed = False

    def refresh_filenames(self):
        self.obsdir = self.ui.lineEdit_obs.text()
        self.fstart = self.ui.fstart.text()
        self.fend = self.ui.fend.text()

    def init_last_selection(self):
        self.ui.lineEdit_obs.setText(self.obsdir)
        self.ui.fstart.setText(self.fstart)
        self.ui.fend.setText(self.fend)

    def get_survey(self):
        return self.survey

    def connectButtons(self):
        QtCore.QObject.connect(self.ui.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseObsdir)

    def chooseObsdir(self):
        self.ui.lineEdit_obs.setText(browseDir('Choose observation directory.'))


class Gen_Survey_from_SR(object):
    def __init__(self, mainwindow):
        self.mainwindow = mainwindow
        self.survey = None
        self.obsdir = None
        self.srcfile = None
        self.recfile = None
        self.fstart = 'shot'
        self.fend = '.dat'
        self.init_dialog()
        self.start_dialog()

    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_generate_survey()
        ui.setupUi(qdialog)
        self.ui = ui
        self.qdialog = qdialog
        self.connectButtons()

    def start_dialog(self):
        self.init_last_selection()
        if self.qdialog.exec_():
            self.refresh_filenames()
            self.survey = activeSeismoPick.Survey(self.obsdir, self.srcfile, self.recfile,
                                                  useDefaultParas = True,
                                                  fstart = self.fstart, fend = self.fend)
            self.executed = True
        else:
            self.refresh_filenames()
            self.executed = False

    def refresh_filenames(self):
        self.obsdir = self.ui.lineEdit_obs.text()
        self.srcfile = self.ui.lineEdit_src.text()
        self.recfile = self.ui.lineEdit_rec.text()
        self.fstart = self.ui.fstart.text()
        self.fend = self.ui.fend.text()

    def init_last_selection(self):
        self.ui.lineEdit_obs.setText(self.obsdir)
        self.ui.lineEdit_src.setText(self.srcfile)
        self.ui.lineEdit_rec.setText(self.recfile)
        self.ui.fstart.setText(self.fstart)
        self.ui.fend.setText(self.fend)

    def get_survey(self):
        return self.survey

    def connectButtons(self):
        QtCore.QObject.connect(self.ui.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseObsdir)
        QtCore.QObject.connect(self.ui.pushButton_src, QtCore.SIGNAL("clicked()"), self.chooseSourcefile)
        QtCore.QObject.connect(self.ui.pushButton_rec, QtCore.SIGNAL("clicked()"), self.chooseRecfile)

    def chooseObsdir(self):
        self.ui.lineEdit_obs.setText(browseDir('Choose observation directory.'))

    def chooseSourcefile(self):
        self.ui.lineEdit_src.setText(openFile('Open sourcefile.'))

    def chooseRecfile(self):
        self.ui.lineEdit_rec.setText(openFile('Open receiverfile.'))



class Call_autopicker(object):
    def __init__(self, mainwindow, survey):
        self.mainwindow = mainwindow
        self.survey = survey
        self.maxSRdist = None
        self.init_dialog()
        self.refresh_selection()                                        
        self.start_dialog()
        
    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_picking_parameters()
        ui.setupUi(qdialog)
        ui.ncores.setMaximum(getMaxCPU())
        self.ui = ui
        self.qdialog = qdialog
        self.initDynSNRplot()
        self.connectButtons()

    def getMaxSRdist(self):
        if self.maxSRdist is not None:
            return self.maxSRdist
        else:
            SRdists = []
            for shot in self.survey.data.values():
                for traceID in shot.getTraceIDlist():
                    SRdists.append(shot.getDistance(traceID))
            self.maxSRdist = max(SRdists)
            return self.maxSRdist
        
    def initDynSNRplot(self):
        self.snrFig = Figure()
        self.snrCanvas = FigureCanvas(self.snrFig)
        self.ui.vlayout_plot.addWidget(self.snrCanvas)
        self.snrToolbar = NavigationToolbar(self.snrCanvas, self.mainwindow)
        self.ui.vlayout_plot.addWidget(self.snrToolbar)

    def plotDynSNR(self):
        fig = self.snrFig
        if fig.axes == []:
            ax = fig.add_subplot(111)
            xlim = None
            ylim = None
        else:
            ax = fig.axes[0]
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
        ax.clear()
        snrthresholds = []
        dists_p = []; snr_p = []
        shiftSNR = float(self.ui.shift_snr.value())
        shiftDist = float(self.ui.shift_dist.value())
        p1 = float(self.ui.p1.value())
        p2 = float(self.ui.p2.value())
        dists = np.arange(0, self.getMaxSRdist() + 1, 1)
        if self.survey.picked:
            for shot in self.survey.data.values():
                for traceID in shot.getTraceIDlist():
                    dists_p.append(shot.getDistance(traceID))
                    snr_p.append(shot.getSNR(traceID)[0])

            ax.scatter(dists_p, snr_p, s = 0.5, c='k')
            
        for dist in dists:
            dist += shiftDist
            snrthresholds.append(surveyUtils.snr_fit_func(surveyUtils.get_fit_fn(p1, p2), dist, shiftSNR))
        ax.plot(dists, snrthresholds)
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('SNR')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        self.snrCanvas.draw()

    def start_dialog(self):
        self.init_last_selection()
        self.plotDynSNR()
        if self.qdialog.exec_():
            self.refresh_selection()

            if self.AIC == True:
                HosAic = 'aic'
            else:
                HosAic = 'hos'

            surveyUtils.setDynamicFittedSNR(self.survey.getShotDict(), shiftdist = self.shiftDist,
                                            shiftSNR = self.shiftSNR, p1 = self.p1, p2 = self.p2)

            self.survey.pickAllShots(vmin = self.vmin, vmax = self.vmax,
                                     folm = self.folm/100., HosAic = HosAic,
                                     aicwindow = self.aicwindow, cores = self.ncores)
            
            #QtGui.qApp.processEvents() # test
            self.executed = True
        else:
            self.refresh_selection()
            self.executed = False

    def refreshFolm(self):
        self.ui.label_folm.setText('%s %%'%self.ui.slider_folm.value())

    def refresh_selection(self):
        self.ncores = int(self.ui.ncores.value())
        self.vmin = float(self.ui.lineEdit_vmin.text())
        self.vmax = float(self.ui.lineEdit_vmax.text())
        self.folm = float(self.ui.slider_folm.value())
        self.AIC = self.ui.checkBox_AIC.isChecked()
        self.aicwindow = (int(self.ui.lineEdit_aicleft.text()), int(self.ui.lineEdit_aicright.text()))
        self.shiftSNR = float(self.ui.shift_snr.value())
        self.shiftDist = float(self.ui.shift_dist.value())
        self.p1 = float(self.ui.p1.value())
        self.p2 = float(self.ui.p2.value())

    def init_last_selection(self):
        self.ui.ncores.setValue(self.ncores)
        self.ui.lineEdit_vmin.setText(str(self.vmin))
        self.ui.lineEdit_vmax.setText(str(self.vmax))
        self.ui.slider_folm.setValue(self.folm)
        self.ui.checkBox_AIC.setChecked(self.AIC)
        self.ui.lineEdit_aicleft.setText(str(self.aicwindow[0]))
        self.ui.lineEdit_aicright.setText(str(self.aicwindow[1]))
        self.ui.shift_snr.setValue(self.shiftSNR)
        self.ui.shift_snr.setValue(self.shiftDist)
        self.ui.p1.setValue(self.p1)
        self.ui.p2.setValue(self.p2)

    def get_survey(self):
        return self.survey

    def connectButtons(self):
        QtCore.QObject.connect(self.ui.slider_folm, QtCore.SIGNAL("valueChanged(int)"), self.refreshFolm)
        QtCore.QObject.connect(self.ui.shift_snr, QtCore.SIGNAL("valueChanged(int)"), self.plotDynSNR)
        QtCore.QObject.connect(self.ui.shift_dist, QtCore.SIGNAL("valueChanged(int)"), self.plotDynSNR)
        QtCore.QObject.connect(self.ui.p1, QtCore.SIGNAL("valueChanged(double)"), self.plotDynSNR)
        QtCore.QObject.connect(self.ui.p2, QtCore.SIGNAL("valueChanged(double)"), self.plotDynSNR)

    def chooseObsdir(self):
        self.ui.lineEdit_obs.setText(browseDir('Choose observation directory.'))

    def chooseSourcefile(self):
        self.ui.lineEdit_src.setText(openFile('Open sourcefile.'))

    def chooseRecfile(self):
        self.ui.lineEdit_rec.setText(openFile('Open receiverfile.'))

