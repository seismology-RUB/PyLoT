#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from PySide import QtCore, QtGui
from pylot.core.active import surveyUtils, activeSeismoPick, seismicArrayPreparation, fmtomoUtils
from generate_survey_layout import Ui_generate_survey
from generate_survey_layout_minimal import Ui_generate_survey_minimal
from generate_seisarray_layout import Ui_generate_seisarray
from picking_parameters_layout import Ui_picking_parameters
from fmtomo_parameters_layout import Ui_fmtomo_parameters
from vtk_tools_layout import Ui_vtk_tools
from postprocessing_layout import Ui_postprocessing

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
        if self.qdialog.exec_():
            self.refresh_selection()
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
            self.refresh_selection()
            self.executed = False

    def refresh_selection(self):
        self.srcfile = self.ui.lineEdit_src.text()
        self.recfile = self.ui.lineEdit_rec.text()
        self.ptsfile = self.ui.lineEdit_pts.text()

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
        if self.qdialog.exec_():
            self.refresh_selection()
            self.survey = activeSeismoPick.Survey(self.obsdir, seisArray = self.seisarray,
                                                  useDefaultParas = True, fstart = self.fstart,
                                                  fend = self.fend)
            self.executed = True
        else:
            self.refresh_selection()
            self.executed = False

    def update_seisarray(self, seisarray):
        self.seisarray = seisarray
        
    def refresh_selection(self):
        self.obsdir = self.ui.lineEdit_obs.text()
        self.fstart = self.ui.fstart.text()
        self.fend = self.ui.fend.text()

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
        if self.qdialog.exec_():
            self.refresh_selection()
            self.survey = activeSeismoPick.Survey(self.obsdir, self.srcfile, self.recfile,
                                                  useDefaultParas = True,
                                                  fstart = self.fstart, fend = self.fend)
            self.executed = True
        else:
            self.refresh_selection()
            self.executed = False

    def refresh_selection(self):
        self.obsdir = self.ui.lineEdit_obs.text()
        self.srcfile = self.ui.lineEdit_src.text()
        self.recfile = self.ui.lineEdit_rec.text()
        self.fstart = self.ui.fstart.text()
        self.fend = self.ui.fend.text()

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
        self.dists_p = []
        self.snr_p = []
        self.lines = []
        self.init_dialog()
        self.refresh_selection()
        self.enableDynSNR(False)                
        self.start_dialog()
        
    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_picking_parameters()
        ui.setupUi(qdialog)
        ui.ncores.setMaximum(getMaxCPU())
        self.ui = ui
        self.qdialog = qdialog
        self.initSNRplot()
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

    def update_survey(self, survey):
        self.survey = survey
        
    def initSNRplot(self):
        self.snrFig = Figure()
        self.snrCanvas = FigureCanvas(self.snrFig)
        self.ui.vlayout_plot.addWidget(self.snrCanvas)
        self.snrToolbar = NavigationToolbar(self.snrCanvas, self.mainwindow)
        self.ui.vlayout_plot.addWidget(self.snrToolbar)

    def prepFigure(self, refresh = True):
        fig = self.snrFig
        if fig.axes == []:
            ax = fig.add_subplot(111)
            xlim = None
            ylim = None
        else:
            ax = fig.axes[0]
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
        #ax.clear()
        if not refresh:
            self.plotPicks(ax)
        else:
            self.clear_lines()

        return fig, ax, xlim, ylim

    def clear_lines(self):
        for line in self.lines:
            line.remove()
        self.lines = []
        
    def finishFigure(self, ax, xlim, ylim):
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    def finishNewFigure(self, ax):
        xlim = None
        ylim = None
        ax.set_xlabel('Distance [m]')
        ax.set_ylabel('SNR')
        return xlim, ylim
        
    def plotPicks(self, ax):
        if self.survey.picked:
            if self.dists_p == [] or self.snr_p == []:
                for shot in self.survey.data.values():
                    for traceID in shot.getTraceIDlist():
                        self.dists_p.append(shot.getDistance(traceID))
                        self.snr_p.append(shot.getSNR(traceID)[0])

            ax.scatter(self.dists_p, self.snr_p, s = 0.5, c='k')
        
    def plotConstSNR(self, refresh = True):
        fig, ax, xlim, ylim = self.prepFigure(refresh)
            
        snrthreshold = float(self.ui.doubleSpinBox_constSNR.text())
        line = ax.hlines(snrthreshold, 0, self.getMaxSRdist(), 'b')
        self.lines.append(line)

        if refresh == False:
            xlim, ylim = self.finishNewFigure(ax)

        if self.survey.picked:
            self.finishFigure(ax, xlim, ylim)
            
        self.snrCanvas.draw()
            
    def plotDynSNR(self, refresh = True):
        fig, ax, xlim, ylim = self.prepFigure(refresh)
        snrthresholds = []
        shiftSNR = float(self.ui.shift_snr.value())
        shiftDist = float(self.ui.shift_dist.value())
        p1 = float(self.ui.p1.value())
        p2 = float(self.ui.p2.value())
        dists = np.arange(0, self.getMaxSRdist() + 1, 1)
        
        for dist in dists:
            dist += shiftDist
            snrthresholds.append(surveyUtils.snr_fit_func(surveyUtils.get_fit_fn(p1, p2), dist, shiftSNR))
        self.lines = ax.plot(dists, snrthresholds, 'b')

        if refresh == False:
            xlim, ylim = self.finishNewFigure(ax)

        if self.survey.picked:
            self.finishFigure(ax, xlim, ylim)
            
        self.snrCanvas.draw()

    def plotSNR(self, refresh = True):
        if self.ui.radioButton_const.isChecked():
            self.plotConstSNR(refresh)
        if self.ui.radioButton_dyn.isChecked():
            self.plotDynSNR(refresh)

    def snr_toggle(self):
        if self.ui.radioButton_const.isChecked():
            self.enableDynSNR(False)
            self.enableConstSNR(True)
        if self.ui.radioButton_dyn.isChecked():
            self.enableConstSNR(False)
            self.enableDynSNR(True)
        self.plotSNR(refresh = True)

    def enableDynSNR(self, bool):
        self.ui.shift_dist.setEnabled(bool)
        self.ui.shift_snr.setEnabled(bool)
        self.ui.p1.setEnabled(bool)
        self.ui.p2.setEnabled(bool)
            
    def enableConstSNR(self, bool):
        self.ui.doubleSpinBox_constSNR.setEnabled(bool)
        
    def start_dialog(self):
        self.plotSNR(refresh = False)
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
            self.clear_lines()
        else:
            self.refresh_selection()
            self.executed = False
            self.clear_lines()

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

    def connectButtons(self):
        QtCore.QObject.connect(self.ui.slider_folm, QtCore.SIGNAL("valueChanged(int)"), self.refreshFolm)
        QtCore.QObject.connect(self.ui.shift_snr, QtCore.SIGNAL("valueChanged(int)"), self.plotSNR)
        QtCore.QObject.connect(self.ui.shift_dist, QtCore.SIGNAL("valueChanged(int)"), self.plotSNR)
        QtCore.QObject.connect(self.ui.p1, QtCore.SIGNAL("valueChanged(double)"), self.plotSNR)
        QtCore.QObject.connect(self.ui.p2, QtCore.SIGNAL("valueChanged(double)"), self.plotSNR)
        QtCore.QObject.connect(self.ui.doubleSpinBox_constSNR, QtCore.SIGNAL("valueChanged(double)"), self.plotSNR)
        QtCore.QObject.connect(self.ui.radioButton_const, QtCore.SIGNAL("clicked()"), self.snr_toggle)
        QtCore.QObject.connect(self.ui.radioButton_dyn, QtCore.SIGNAL("clicked()"), self.snr_toggle)

    def chooseObsdir(self):
        self.ui.lineEdit_obs.setText(browseDir('Choose observation directory.'))

    def chooseSourcefile(self):
        self.ui.lineEdit_src.setText(openFile('Open sourcefile.'))

    def chooseRecfile(self):
        self.ui.lineEdit_rec.setText(openFile('Open receiverfile.'))


class Call_FMTOMO(object):
    def __init__(self, mainwindow, survey):
        self.mainwindow = mainwindow
        self.survey = survey
        self.init_dialog()
        self.refresh_selection()                                        
        self.start_dialog()

    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_fmtomo_parameters()
        ui.setupUi(qdialog)
        ui.nproc.setMaximum(getMaxCPU())
        self.ui = ui
        self.qdialog = qdialog
        self.connectButtons()
        
    def start_dialog(self):
        if self.qdialog.exec_():
            self.refresh_selection()

            if not os.path.isdir(self.picks_dir):
                err = os.mkdir(self.picks_dir) # error not handled yet

            self.survey.exportFMTOMO(self.picks_dir)

            cwd = os.getcwd()
            interpolationMethod = 'linear'
            os.chdir(self.simuldir)
            if self.survey.seisarray.twoDim:
                interpolationMethod = 'nearest'
            self.survey.seisarray.generateFMTOMOinputFromArray(self.propgrid, self.vgrid, (self.bbot, self.btop),
                                                               self.cushionfactor/100., interpolationMethod,
                                                               customgrid = self.customgrid, writeVTK = True)
            os.chdir(cwd)

            tomo = fmtomoUtils.Tomo3d(self.fmtomo_dir, self.simuldir)
            tomo.runTOMO3D(self.nproc, self.nIter)

            self.executed = True
        else:
            self.refresh_selection()
            self.executed = False

    def update_survey(self, survey):
        self.survey = survey
        
    def refresh_selection(self):
        self.fmtomo_dir = self.ui.fmtomo_dir.text()
        self.nIter = int(self.ui.nIter.value())
        self.nproc = int(self.ui.nproc.value())
        self.btop = float(self.ui.btop.text())
        self.bbot = float(self.ui.bbot.text())
        self.propgrid = (int(self.ui.pgrid_x.value()), int(self.ui.pgrid_y.value()), int(self.ui.pgrid_z.value()))
        self.vgrid = (int(self.ui.invgrid_x.value()), int(self.ui.invgrid_y.value()), int(self.ui.invgrid_z.value()))
        self.cushionfactor = float(self.ui.cushion.value())
        self.customgrid = self.ui.customgrid.text()
        self.simuldir = self.ui.simuldir.text()
        self.picks_dir = os.path.join(self.simuldir, 'picks')
        
    def connectButtons(self):
        QtCore.QObject.connect(self.ui.browse_tomodir, QtCore.SIGNAL("clicked()"), self.chooseFMTOMOdir)
        QtCore.QObject.connect(self.ui.browse_customgrid, QtCore.SIGNAL("clicked()"), self.chooseCustomgrid)
        QtCore.QObject.connect(self.ui.browse_simuldir, QtCore.SIGNAL("clicked()"), self.chooseSimuldir)

    def chooseFMTOMOdir(self):
        self.ui.fmtomo_dir.setText(browseDir())

    def chooseCustomgrid(self):
        self.ui.customgrid.setText(openFile())

    def chooseSimuldir(self):
        self.ui.simuldir.setText(browseDir())


class Call_VTK_dialog(object):
    def __init__(self, mainwindow):
        self.mainwindow = mainwindow
        self.init_dialog()
        self.refresh_selection()                                        
        self.start_dialog()

    def init_dialog(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_vtk_tools()
        ui.setupUi(qdialog)
        self.ui = ui
        self.qdialog = qdialog
        self.connectButtons()
        
    def start_dialog(self):
        self.qdialog.exec_()
        self.refresh_selection()

    def refresh_selection(self):
        self.vg = self.ui.lineEdit_vg.text()
        self.vgout = self.ui.lineEdit_vgout.text()
        self.rays = self.ui.lineEdit_rays.text()
        self.raysout = self.ui.lineEdit_raysout.text()
        
    def checkVgStartButton(self):
        ui = self.ui
        if ui.radioButton_rel.isChecked():
            if ui.lineEdit_vg.text() != '' and ui.lineEdit_vgref.text() != '':
                ui.start_vg.setEnabled(True)
            else:
                ui.start_vg.setEnabled(False)
        if ui.radioButton_abs.isChecked():
            if ui.lineEdit_vg.text() != '':
                ui.start_vg.setEnabled(True)
            else:
                ui.start_vg.setEnabled(False)                

    def checkRaysStartButton(self):
        ui = self.ui
        if ui.lineEdit_rays.text() != '' and ui.lineEdit_raysout.text() != '':
            ui.start_rays.setEnabled(True)
        else:
            ui.start_rays.setEnabled(False)

    def connectButtons(self):
        QtCore.QObject.connect(self.ui.pushButton_vg, QtCore.SIGNAL("clicked()"), self.chooseVgrid)
        QtCore.QObject.connect(self.ui.pushButton_vgref, QtCore.SIGNAL("clicked()"), self.chooseVgridref)
        QtCore.QObject.connect(self.ui.pushButton_rays, QtCore.SIGNAL("clicked()"), self.chooseRaysIn)
        QtCore.QObject.connect(self.ui.pushButton_raysout, QtCore.SIGNAL("clicked()"), self.chooseRaysOutDir)
        QtCore.QObject.connect(self.ui.pushButton_vtkout, QtCore.SIGNAL("clicked()"), self.newFileVTK)
        QtCore.QObject.connect(self.ui.pushButton_parav, QtCore.SIGNAL("clicked()"), self.openFileParaview)
        QtCore.QObject.connect(self.ui.start_vg, QtCore.SIGNAL("clicked()"), self.startvgvtk)
        QtCore.QObject.connect(self.ui.start_rays, QtCore.SIGNAL("clicked()"), self.startraysvtk)
        QtCore.QObject.connect(self.ui.radioButton_rel, QtCore.SIGNAL("clicked()"), self.activateVgref)
        QtCore.QObject.connect(self.ui.radioButton_abs, QtCore.SIGNAL("clicked()"), self.deactivateVgref)

    def openFileParaview(self):
        os.system('paraview %s &'%self.ui.lineEdit_vgout.text())

    def activateVgref(self):
        self.ui.lineEdit_vgref.setEnabled(True)
        self.ui.pushButton_vgref.setEnabled(True)

    def deactivateVgref(self):
        self.ui.lineEdit_vgref.setEnabled(False)
        self.ui.pushButton_vgref.setEnabled(False)

    def chooseVgrid(self):
        self.ui.lineEdit_vg.setText(openFile())
        self.checkVgStartButton()

    def chooseVgridref(self):
        self.ui.lineEdit_vgref.setText(openFile())
        self.checkVgStartButton()

    def chooseRaysIn(self):
        self.ui.lineEdit_rays.setText(openFile())
        self.checkRaysStartButton()

    def chooseRaysOutDir(self):
        self.ui.lineEdit_raysout.setText(browseDir())
        self.checkRaysStartButton()

    def startvgvtk(self):
        ui = self.ui
        if ui.lineEdit_vgout.text() == '':
            return            
        if ui.radioButton_abs.isChecked():
            fmtomoUtils.vgrids2VTK(inputfile = ui.lineEdit_vg.text(),
                                   outputfile = ui.lineEdit_vgout.text(),
                                   absOrRel='abs')
        elif ui.radioButton_rel.isChecked():
            fmtomoUtils.vgrids2VTK(inputfile = ui.lineEdit_vg.text(),
                                   outputfile = ui.lineEdit_vgout.text(),
                                   absOrRel='rel',
                                   inputfileref = ui.lineEdit_vgref.text())

    def startraysvtk(self):
        ui = self.ui
        fmtomoUtils.rays2VTK(ui.lineEdit_rays.text(), ui.lineEdit_raysout.text())

    def newFileVTK(self):
        self.ui.lineEdit_vgout.setText(saveFile())

class Postprocessing(object):
    def __init__(self, mainwindow, survey):
        self.mainwindow = mainwindow
        self.survey = survey
        self.init_dialog()
        self.start_dialog()

    def init_dialog(self):
        qwidget = QtGui.QWidget()#
        ui = Ui_postprocessing()
        ui.setupUi(qwidget)
        self.ui = ui
        self.qwidget = qwidget
        self.initPlot()
        self.plot()
        #self.connectButtons()

    def start_dialog(self):
        self.qwidget.show()
        # if self.qwidget.exec_():
        #     #self.refresh_selection()
        #     self.executed = True
        # else:
        #     self.refresh_selection()
        #     self.executed = False

    def initPlot(self):
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.ui.verticalLayout_plot.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self.mainwindow)
        self.ui.verticalLayout_plot.addWidget(self.toolbar)
        
    def plot(self):
        survey = self.survey
        ax = self.figure.add_subplot(111)
        dist, pick, snrlog, pickerror, spe = survey.preparePlotAllPicks(plotRemoved = False)
        ax, cbar, sc = survey.createPlot(dist, pick, snrlog, '123', ax = ax, cbar = None)
        self.cbar = self.figure.colorbar(sc, fraction=0.05)
        self.ax = ax
        
    # def refresh_selection(self):
    #     self.obsdir = self.ui.lineEdit_obs.text()
    #     self.fstart = self.ui.fstart.text()
    #     self.fend = self.ui.fend.text()

    def update_survey(self, survey):
        self.survey = survey
        
    def get_survey(self):
        return self.survey

    # def connectButtons(self):
    #     QtCore.QObject.connect(self.ui.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseObsdir)

    # def chooseObsdir(self):
    #     self.ui.lineEdit_obs.setText(browseDir('Choose observation directory.'))
