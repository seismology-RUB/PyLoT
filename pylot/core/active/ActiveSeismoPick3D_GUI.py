#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'

from PySide import QtCore, QtGui, QtCore
from asp3d_layout import *
from fmtomo_parameters_layout import *
from generate_survey_layout import *
from generate_survey_layout_minimal import *
from generate_seisarray_layout import *
from picking_parameters_layout import *
from pylot.core.active import activeSeismoPick, surveyUtils, fmtomoUtils, seismicArrayPreparation

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

class gui_control(object):
    def __init__(self):
        self.mainwindow = MainWindow
        self.mainUI = ui
        self.connectButtons()
        self.survey = None
        self.seisarray = None
        self.seisArrayFigure = None
        self.cancelpixmap = self.mainwindow.style().standardPixmap(QtGui.QStyle.SP_DialogCancelButton)
        self.applypixmap = self.mainwindow.style().standardPixmap(QtGui.QStyle.SP_DialogApplyButton)
        self.setInitStates()
        self.addArrayPlot()
        self.addStatPlots()

    def setInitStates(self):
        self.setPickState(False)
        self.setSurveyState(False)
        self.setSeisArrayState(False)
        self.setConnected2SurveyState(False)

    def connectButtons(self):
        QtCore.QObject.connect(self.mainUI.gen_new_seisarray, QtCore.SIGNAL("clicked()"), self.gen_seisarray)
        QtCore.QObject.connect(self.mainUI.load_seisarray, QtCore.SIGNAL("clicked()"), self.load_seisarray)
        QtCore.QObject.connect(self.mainUI.save_seisarray, QtCore.SIGNAL("clicked()"), self.save_seisarray)
        QtCore.QObject.connect(self.mainUI.connect_to_survey, QtCore.SIGNAL("clicked()"), self.connect2Survey)
        QtCore.QObject.connect(self.mainUI.interpolate_receivers, QtCore.SIGNAL("clicked()"), self.interpolate_receivers)
        QtCore.QObject.connect(self.mainUI.gen_new_survey, QtCore.SIGNAL("clicked()"), self.gen_survey)
        QtCore.QObject.connect(self.mainUI.load_survey, QtCore.SIGNAL("clicked()"), self.load_survey)
        QtCore.QObject.connect(self.mainUI.save_survey, QtCore.SIGNAL("clicked()"), self.save_survey)
        QtCore.QObject.connect(self.mainUI.picker, QtCore.SIGNAL("clicked()"), self.callPicker)
        QtCore.QObject.connect(self.mainUI.postprocessing, QtCore.SIGNAL("clicked()"), self.postprocessing)
        QtCore.QObject.connect(self.mainUI.fmtomo, QtCore.SIGNAL("clicked()"), self.startFMTOMO)
        QtCore.QObject.connect(self.mainUI.comboBox, QtCore.SIGNAL("activated(int)"), self.replotStat)

    def gen_seisarray(self):
        disconnect = False
        if self.checkSeisArrayState():
            if not self.continueDialogExists('Seismic Array'):
                return
        if self.checkConnected2SurveyState():
            if not self.continueDialogMessage('Seismic Array connected to present Survey.\n'
                                              'Continuation will disconnect the Seismic Array.'):
                return
            else:
                self.survey.seisarray = None
                disconnect = True
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_generate_seisarray()
        ui.setupUi(qdialog)
        self.gen_new_seisarray = ui
        self.connectButtons_gen_seisarray()
        if qdialog.exec_():
            srcfile = self.gen_new_seisarray.lineEdit_src.text()
            recfile = self.gen_new_seisarray.lineEdit_rec.text()
            ptsfile = self.gen_new_seisarray.lineEdit_pts.text()
            self.seisarray = seismicArrayPreparation.SeisArray(recfile)
            if len(srcfile) > 0:
                self.seisarray.addSourceLocations(srcfile)
            if len(ptsfile) > 0:
                self.seisarray.addMeasuredTopographyPoints(ptsfile)
            if disconnect:
                self.setConnected2SurveyState(False)
            self.setSeisArrayState(True)

    def gen_survey(self):
        if self.checkSurveyState():
            if not self.continueDialogExists('Survey'):
                return
        if self.checkSeisArrayState():
            if self.continueDialogMessage('Use geometry information of active Seismic Array?'):
               if self.gen_survey_fromSeisArray():
                   self.initNewSurvey()
                   return
               else:
                   return
        if self.gen_survey_fromSRfiles():
            self.initNewSurvey()

    def initNewSurvey(self):
        self.survey.setArtificialPick(0, 0) # artificial pick at source origin                                         
        surveyUtils.setDynamicFittedSNR(self.survey.getShotDict())
        self.setSurveyState(True)
        self.setPickState(False)

    def gen_survey_fromSeisArray(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_generate_survey_minimal()
        ui.setupUi(qdialog)
        self.gen_new_survey_min = ui
        self.connectButtons_gen_survey_min()
        if qdialog.exec_():
            obsdir = self.gen_new_survey_min.lineEdit_obs.text()
            self.survey = activeSeismoPick.Survey(obsdir, seisArray = self.seisarray,
                                                  useDefaultParas = True)
            self.setConnected2SurveyState(True)
            return True

    def gen_survey_fromSRfiles(self):
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_generate_survey()
        ui.setupUi(qdialog)
        self.gen_new_survey = ui
        self.connectButtons_gen_survey()
        if qdialog.exec_():
            srcfile = self.gen_new_survey.lineEdit_src.text()
            recfile = self.gen_new_survey.lineEdit_rec.text()
            obsdir = self.gen_new_survey.lineEdit_obs.text()
            self.survey = activeSeismoPick.Survey(obsdir, srcfile, recfile,
                                                  useDefaultParas = True)
            self.setConnected2SurveyState(False)
            return True

    def addArrayPlot(self):
        self.seisArrayFigure = Figure()
        self.seisArrayCanvas = FigureCanvas(self.seisArrayFigure)
        self.mainUI.horizontalLayout_tr.addWidget(self.seisArrayCanvas)

    def addStatPlots(self):
        self.statFigure_left = Figure()
        self.statCanvas_left = FigureCanvas(self.statFigure_left)
        self.mainUI.horizontalLayout_br.addWidget(self.statCanvas_left)
        self.statFigure_right = Figure()
        self.statCanvas_right = FigureCanvas(self.statFigure_right)
        self.mainUI.horizontalLayout_br.addWidget(self.statCanvas_right)
        self.addItems2ComboBox()

    def addItems2ComboBox(self):
        self.mainUI.comboBox.insertItem(0, 'picked traces')
        self.mainUI.comboBox.insertItem(1, 'mean SNR')
        self.mainUI.comboBox.insertItem(2, 'median SNR')
        self.mainUI.comboBox.insertItem(3, 'mean SPE')
        self.mainUI.comboBox.insertItem(4, 'median SPE')
        self.mainUI.comboBox.setEnabled(False)

    def addArrayAxes(self):
        self.seisArrayAx = self.seisArrayFigure.add_subplot(111)

    def addStatAxes(self):
        self.statAx_left = self.statFigure_left.add_subplot(111)
        self.statAx_right = self.statFigure_right.add_subplot(111)

    def replotArray(self):
        self.seisArrayFigure.clf()
        self.addArrayAxes()
        self.plotArray()
        self.seisArrayCanvas.draw()

    def plotArray(self):
        self.seisarray.plotArray2D(self.seisArrayAx, highlight_measured = True)

    def replotStat(self):
        self.statFigure_left.clf()
        self.statFigure_right.clf()
        self.addStatAxes()
        self.plotStat()
        self.statCanvas_left.draw()
        self.statCanvas_right.draw()

    def plotStat(self):
        if self.checkPickState():
            surveyUtils.plotScatterStats4Receivers(self.survey, self.mainUI.comboBox.currentText(), self.statAx_left)
            surveyUtils.plotScatterStats4Shots(self.survey, self.mainUI.comboBox.currentText(), self.statAx_right)

    def interpolate_receivers(self):
        if not self.checkSeisArrayState():
            self.printDialogMessage('No Seismic Array defined.')
            return
        self.seisarray.interpolateAll()
        self.replotArray()

    def getPickParameters(self, ui, Picking_parameters):
        if Picking_parameters.exec_():
            ncores = int(ui.ncores.value())
            vmin = float(ui.lineEdit_vmin.text())
            vmax = float(ui.lineEdit_vmax.text())
            folm = float(ui.lineEdit_folm.text())
            AIC = ui.checkBox_AIC.isChecked()
            aicwindow = (int(ui.lineEdit_aicleft.text()), int(ui.lineEdit_aicright.text()))
            return ncores, vmin, vmax, folm, AIC, aicwindow

    def connect2Survey(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        if not self.checkSeisArrayState():
            self.printDialogMessage('Got no Seismic Array.')
            return
        if self.checkConnected2SurveyState():
            if not self.continueDialogMessage('Existing Survey already got Seismic Array object. Continue?'):
                return
        self.survey.seisarray = self.seisarray
        self.setConnected2SurveyState(True)
        self.survey._initiate_SRfiles()
        print('Connected Seismic Array to active Survey object.')

    def getMaxCPU(self):
        import multiprocessing
        return multiprocessing.cpu_count()

    def callPicker(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        if self.checkPickState():
            if not self.continueDialogMessage('Survey already picked. Continue?'):
                return
        Picking_parameters = QtGui.QDialog(self.mainwindow)
        ui = Ui_picking_parameters()
        ui.setupUi(Picking_parameters)
        ui.ncores.setMaximum(self.getMaxCPU())
        try:
            ncores, vmin, vmax, folm, AIC, aicwindow = self.getPickParameters(ui, Picking_parameters)
        except TypeError:
            return
        
        if AIC == True:
            HosAic = 'aic'
        else:
            HosAic = 'hos'

        self.survey.pickAllShots(vmin = vmin, vmax = vmax,
                                 folm = folm, HosAic = HosAic,
                                 aicwindow = aicwindow, cores = ncores)
        self.setPickState(True)

    def startFMTOMO(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        fmtomo_parameters = QtGui.QDialog(self.mainwindow)
        ui = Ui_fmtomo_parameters()
        ui.setupUi(fmtomo_parameters)                
        ui.nproc.setMaximum(self.getMaxCPU())

        self.fmtomo_parameters_ui = ui
        self.connectButtons_startFMTOMO()
        self.getFMTOMOparameters(ui, fmtomo_parameters)

    def getFMTOMOparameters(self, ui, fmtomo_parameters):
        if fmtomo_parameters.exec_():
            fmtomo_dir = ui.fmtomo_dir.text()
            nIter = int(ui.nIter.value())
            nproc = int(ui.nproc.value())
            btop = float(ui.btop.text())
            bbot = float(ui.bbot.text())
            propgrid = (int(ui.pgrid_x.text()), int(ui.pgrid_y.text()), int(ui.pgrid_z.text()))
            vgrid = (int(ui.invgrid_x.text()), int(ui.invgrid_y.text()), int(ui.invgrid_z.text()))
            cushionfactor = float(ui.cushion.text())/100.
            customgrid = ui.customgrid.text()
            simuldir = ui.simuldir.text()
            picks_dir = os.path.join(simuldir, 'picks')

            if not os.path.isdir(picks_dir):
                err = os.mkdir(picks_dir)

            self.survey.exportFMTOMO(picks_dir)

            cwd = os.getcwd()
            interpolationMethod = 'linear'
            os.chdir(simuldir)
            self.survey.seisarray.generateFMTOMOinputFromArray(propgrid, vgrid, (bbot, btop), cushionfactor,
                                                          interpolationMethod, customgrid = customgrid,
                                                          writeVTK = False)
            os.chdir(cwd)

            tomo = fmtomoUtils.Tomo3d(fmtomo_dir, simuldir)
            tomo.runTOMO3D(nproc, nIter)
                            
    def connectButtons_startFMTOMO(self):
        QtCore.QObject.connect(self.fmtomo_parameters_ui.browse_tomodir, QtCore.SIGNAL("clicked()"), self.chooseFMTOMOdir)
        QtCore.QObject.connect(self.fmtomo_parameters_ui.browse_customgrid, QtCore.SIGNAL("clicked()"), self.chooseCustomgrid)
        QtCore.QObject.connect(self.fmtomo_parameters_ui.browse_simuldir, QtCore.SIGNAL("clicked()"), self.chooseSimuldir)

    def chooseFMTOMOdir(self):
        self.fmtomo_parameters_ui.fmtomo_dir.setText(self.browseDir())

    def chooseCustomgrid(self):
        self.fmtomo_parameters_ui.customgrid.setText(self.openFile())

    def chooseSimuldir(self):
        self.fmtomo_parameters_ui.simuldir.setText(self.browseDir())

    def chooseSeisarray(self):
        self.fmtomo_parameters_ui.seisarray.setText(self.openFile())

    def postprocessing(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        self.survey.plotAllPicks()
        
    def load_survey(self):
        if self.checkSurveyState():
            if not self.continueDialogExists('Survey'):
                return
        filename = self.openFile()
        if filename is None:
            return
        try:
            survey = activeSeismoPick.Survey.from_pickle(filename)
        except:
            self.printDialogMessage('Could not load object %s.'%filename)
            return
        if not type(survey) == activeSeismoPick.Survey:
            self.printDialogMessage('Wrong input file of type %s, expected %s.'
                  %(type(survey), activeSeismoPick.Survey))
            return
        if self.checkSeisArrayState() and survey.seisarray is not None:
            if not self.continueDialogMessage('Survey got existing Seismic Array.'
                                              ' Do you want to overwrite the current Seismic Array?'):
                return
        self.survey = survey
        self.setSurveyState(True)
        if self.survey.picked:
            self.setPickState(True)
        else:
            self.setPickState(False)
        if self.survey.seisarray != None:
            self.seisarray = self.survey.seisarray
            self.setConnected2SurveyState(True)
            self.setSeisArrayState(True)
            self.printDialogMessage('Loaded Survey with active Seismic Array.')
        else:
            self.setConnected2SurveyState(False)
            self.setSeisArrayState(False)
            self.printDialogMessage('Loaded Survey.')

    def load_seisarray(self):
        disconnect = False
        if self.checkSeisArrayState():
            if not self.continueDialogExists('Seismic Array'):
                return
        if self.checkConnected2SurveyState():
            if not self.continueDialogMessage('Seismic Array connected to present Survey.\n'
                                              'Continuation will disconnect the Seismic Array.'):
                return
            else:
                self.survey.seisarray = None
                disconnect = True

        filename = self.openFile()
        if filename is None:
            return
        try:
            seisarray = seismicArrayPreparation.SeisArray.from_pickle(filename)
        except:
            self.printDialogMessage('Could not load object %s.'%filename)
            return
        if not type(seisarray) == seismicArrayPreparation.SeisArray:
            self.printDialogMessage('Wrong input file of type %s, expected %s.'
                  %(type(survey), seismicArrayPreparation.SeisArray))
            return
        if disconnect:
            self.setConnected2SurveyState(False)
        self.seisarray = seisarray
        self.replotArray()
        self.setSeisArrayState(True)

    def save_seisarray(self):
        if not self.checkSeisArrayState():
            self.printDialogMessage('No Seismic Array defined.')
            return
        filename = self.saveFile()
        if filename is None:
            return
        self.seisarray.saveSeisArray(filename)

    def save_survey(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        filename = self.saveFile()
        if filename is None:
            return
        self.survey.saveSurvey(filename)

    def setSurveyState(self, state):
        if state == True:
            self.mainUI.survey_active.setPixmap(self.applypixmap)
        elif state == False:
            self.mainUI.survey_active.setPixmap(self.cancelpixmap)

    def checkSurveyState(self):
        if self.survey == None:
            return False
        else:
            return True

    def checkSeisArrayState(self):
        if self.seisarray == None:
            return False
        else:
            return True

    def setPickState(self, state):
        if state == True and self.checkSurveyState():
            self.mainUI.picked_active.setPixmap(self.applypixmap)
            self.replotStat()
            self.mainUI.comboBox.setEnabled(True)
            self.survey.picked = True
        elif state == True and self.checkSurveyState() is False:
            self.printDialogMessage('No Survey defined.')
            return
        elif state == False:
            self.mainUI.picked_active.setPixmap(self.cancelpixmap)
            if self.checkSurveyState():
                self.statFigure_left.clf()
                self.statFigure_right.clf()
                self.mainUI.comboBox.setEnabled(False)
                self.survey.picked = False

    def setSeisArrayState(self, state):
        if state == True:
            self.mainUI.seisarray_active.setPixmap(self.applypixmap)
            self.replotArray()
        elif state == False:
            self.mainUI.seisarray_active.setPixmap(self.cancelpixmap)
            if self.seisArrayFigure is not None:
                self.seisArrayFigure.clf()

    def setConnected2SurveyState(self, state):
        if state == True:
            self.mainUI.seisarray_on_survey_active.setPixmap(self.applypixmap)
        elif state == False:
            self.mainUI.seisarray_on_survey_active.setPixmap(self.cancelpixmap)

    def checkConnected2SurveyState(self):
        if self.checkSurveyState():
            if self.survey.seisarray != None:
                return True
        else:
            return False

    def checkPickState(self):
        if not self.survey:
            self.printDialogMessage('No Survey defined.')
            return
        return self.survey.picked

    def printDialogMessage(self, message):
        qmb = QtGui.QMessageBox()
        qmb.setText(message)
        qmb.setStandardButtons(QtGui.QMessageBox.Ok)
        qmb.setIcon(QtGui.QMessageBox.Warning)
        qmb.exec_()

    def continueDialogExists(self, name):
        qmb = QtGui.QMessageBox()
        qmb.setText('%s object already exists. Overwrite?'%name)
        qmb.setStandardButtons(QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
        qmb.setIcon(QtGui.QMessageBox.Warning)
        answer = qmb.exec_()
        if answer == 1024:
            return True
        else:
            return False

    def continueDialogMessage(self, message):
        qmb = QtGui.QMessageBox()
        qmb.setText(message)
        qmb.setStandardButtons(QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
        qmb.setIcon(QtGui.QMessageBox.Warning)
        answer = qmb.exec_()
        if answer == 1024:
            return True
        else:
            return False

    def connectButtons_gen_survey(self):
        QtCore.QObject.connect(self.gen_new_survey.pushButton_rec, QtCore.SIGNAL("clicked()"), self.chooseReceiverfile)
        QtCore.QObject.connect(self.gen_new_survey.pushButton_src, QtCore.SIGNAL("clicked()"), self.chooseSourcefile)
        QtCore.QObject.connect(self.gen_new_survey.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseObsdir)

    def connectButtons_gen_survey_min(self):
        QtCore.QObject.connect(self.gen_new_survey_min.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseObsdir_min)

    def connectButtons_gen_seisarray(self):
        QtCore.QObject.connect(self.gen_new_seisarray.pushButton_rec, QtCore.SIGNAL("clicked()"), self.chooseMeasuredRec)
        QtCore.QObject.connect(self.gen_new_seisarray.pushButton_src, QtCore.SIGNAL("clicked()"), self.chooseMeasuredSrc)
        QtCore.QObject.connect(self.gen_new_seisarray.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseMeasuredPts)

    def chooseMeasuredSrc(self):
        self.gen_new_seisarray.lineEdit_src.setText(self.openFile('Open measured sources file.'))

    def chooseMeasuredRec(self):
        self.gen_new_seisarray.lineEdit_rec.setText(self.openFile('Open measured receivers file.'))

    def chooseMeasuredPts(self):
        self.gen_new_seisarray.lineEdit_pts.setText(self.openFile('Open measured points file.'))

    def chooseSourcefile(self):
        self.gen_new_survey.lineEdit_src.setText(self.openFile('Open sourcefile.'))

    def chooseReceiverfile(self):
        self.gen_new_survey.lineEdit_rec.setText(self.openFile('Open receiverfile.'))

    def chooseObsdir(self):
        self.gen_new_survey.lineEdit_obs.setText(self.browseDir('Choose observation directory.'))

    def chooseObsdir_min(self):
        self.gen_new_survey_min.lineEdit_obs.setText(self.browseDir('Choose observation directory.'))

    def openFile(self, name = 'Open'):
        dialog = QtGui.QFileDialog()
        dialog.setWindowTitle(name)                #not working yet
        filename = dialog.getOpenFileName()
        if len(filename[0]) > 0:
            return filename[0]

    def saveFile(self, name = 'Save'):
        dialog = QtGui.QFileDialog()
        dialog.setWindowTitle(name)
        filename = dialog.getSaveFileName()
        if len(filename[0]) > 0:
            return filename[0]

    def browseDir(self, name = 'Open Directory'):
        dialog = QtGui.QFileDialog()
        dialog.setWindowTitle(name)
        directory = dialog.getExistingDirectory()
        if len(directory) > 0:
            return directory

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    gui = gui_control()
    sys.exit(app.exec_())


