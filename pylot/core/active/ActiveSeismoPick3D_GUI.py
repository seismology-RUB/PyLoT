#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'

from PySide import QtCore, QtGui
from pylot.core.active import activeSeismoPick, surveyUtils, fmtomoUtils, seismicArrayPreparation
from pylot.core.active.gui.asp3d_layout import *
from pylot.core.active.gui.windows import Gen_SeisArray, Gen_Survey_from_SA, Gen_Survey_from_SR, Call_autopicker, Call_FMTOMO, Call_VTK_dialog, Postprocessing
from pylot.core.active.gui.windows import openFile, saveFile, browseDir, getMaxCPU

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
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
        self.addArrayPlot()
        self.addSurfacePlot()
        self.addStatPlots()
        self.setInitStates()
        self.mainUI.progressBar.setVisible(False)
        self.printSurveyTextbox()
        self.printSeisArrayTextbox()
        self.initWindowObjects()

    def initWindowObjects(self):
        self.gsa = None
        self.gssa = None
        self.gssr = None
        self.autopicker = None
        self.fmtomo = None
        self.vtktools = None
        self.postprocessing = None

    def setInitStates(self):
        self.setPickState(False)
        self.setSurveyState(False)
        self.setSeisArrayState(False)
        self.setConnected2SurveyState(False)

    def connectButtons(self):
        QtCore.QObject.connect(self.mainUI.actionGenerate_new_Seismic_Array, QtCore.SIGNAL("triggered()"), self.gen_seisarray)
        QtCore.QObject.connect(self.mainUI.actionLoad_Seismic_Array, QtCore.SIGNAL("triggered()"), self.load_seisarray)
        QtCore.QObject.connect(self.mainUI.actionSave_Seismic_Array, QtCore.SIGNAL("triggered()"), self.save_seisarray)
        QtCore.QObject.connect(self.mainUI.actionLoad_Survey, QtCore.SIGNAL("triggered()"), self.load_survey)
        QtCore.QObject.connect(self.mainUI.actionSave_Survey, QtCore.SIGNAL("triggered()"), self.save_survey)
        QtCore.QObject.connect(self.mainUI.actionConnect_to_Survey, QtCore.SIGNAL("triggered()"), self.connect2Survey)
        QtCore.QObject.connect(self.mainUI.actionInterpolate_Receivers, QtCore.SIGNAL("triggered()"), self.interpolate_receivers)
        QtCore.QObject.connect(self.mainUI.actionGenerate_new_Survey, QtCore.SIGNAL("triggered()"), self.gen_survey)
        QtCore.QObject.connect(self.mainUI.actionAutomatic_Picking, QtCore.SIGNAL("triggered()"), self.startPicker)
        QtCore.QObject.connect(self.mainUI.actionPostprocessing, QtCore.SIGNAL("triggered()"), self.postprocessing)
        QtCore.QObject.connect(self.mainUI.actionStart_FMTOMO_Simulation, QtCore.SIGNAL("triggered()"), self.startFMTOMO)
        QtCore.QObject.connect(self.mainUI.actionVTK_Visualization, QtCore.SIGNAL("triggered()"), self.startVTKtools)
        QtCore.QObject.connect(self.mainUI.actionExit, QtCore.SIGNAL("triggered()"), self.exitApp)
        QtCore.QObject.connect(self.mainUI.actionFullscreen, QtCore.SIGNAL("triggered()"), self.fullscreen)
        QtCore.QObject.connect(self.mainUI.comboBox_stats, QtCore.SIGNAL("activated(int)"), self.refreshPickedWidgets)
        QtCore.QObject.connect(self.mainUI.shot_left, QtCore.SIGNAL("clicked()"), self.decreaseShotnumber)
        QtCore.QObject.connect(self.mainUI.shot_right, QtCore.SIGNAL("clicked()"), self.increaseShotnumber)
        QtCore.QObject.connect(self.mainUI.plot_shot, QtCore.SIGNAL("clicked()"), self.plotShot)

    def fullscreen(self):
        if self.mainUI.actionFullscreen.isChecked():
            MainWindow.showFullScreen()
        else:
            MainWindow.showNormal()                    
        
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
        
        if self.gsa is None:
            self.gsa = Gen_SeisArray(self.mainwindow)
        else:
            self.gsa.start_dialog()

        if self.gsa.executed:
            self.seisarray = self.gsa.get_seisarray()
            if disconnect:
                self.setConnected2SurveyState(False)
            self.setSeisArrayState(True)

    def gen_survey(self):
        if self.checkSurveyState():
            if not self.continueDialogExists('Survey'):
                return
        if self.checkSeisArrayState():
            if len(self.seisarray.getSourceCoordinates()) > 0:
                if self.continueDialogMessage('Use geometry information of active Seismic Array?'):
                    if self.gssa is None:
                        self.gssa = Gen_Survey_from_SA(self.mainwindow, self.seisarray)
                    else:
                        self.gssa.start_dialog()
                        self.update_seisarray(self.seisarray)
                    if self.gssa.executed:
                        self.survey = self.gssa.get_survey()
                        self.initNewSurvey()
                        self.setConnected2SurveyState(True)
                        self.setPickState(False)
                    return
            else:
                if not self.continueDialogMessage('Can not use current Seismic Array,'
                                                  ' because there are no sources given.'):
                    return
        if self.gssr is None:
            self.gssr = Gen_Survey_from_SR(self.mainwindow)
        else:
            self.gssr.start_dialog()
        if self.gssr.executed:
            self.survey = self.gssr.get_survey()
            self.seisarray = self.survey.seisarray
            self.initNewSurvey()
            self.setSeisArrayState(True)
            self.setConnected2SurveyState(True)

    def initNewSurvey(self):
        self.survey.setArtificialPick(0, 0) # artificial pick at source origin
        self.setSurveyState(True)
        self.setPickState(False)

    def addArrayPlot(self):
        self.seisArrayFigure = Figure()
        self.seisArrayCanvas = FigureCanvas(self.seisArrayFigure)
        self.mainUI.verticalLayout_tr1.addWidget(self.seisArrayCanvas)
        self.seisArrayToolbar = NavigationToolbar(self.seisArrayCanvas, self.mainwindow)
        self.mainUI.verticalLayout_tr1.addWidget(self.seisArrayToolbar)

    def addSurfacePlot(self):
        self.surfaceFigure = Figure()
        self.surfaceCanvas = FigureCanvas(self.surfaceFigure)
        self.mainUI.horizontalLayout_tr.addWidget(self.surfaceCanvas)

    def addStatPlots(self):
        self.statFigure_left = Figure()
        self.statCanvas_left = FigureCanvas(self.statFigure_left)
        self.mainUI.verticalLayout_br1.addWidget(self.statCanvas_left)
        self.statToolbar_left = NavigationToolbar(self.statCanvas_left, self.mainwindow)
        self.mainUI.verticalLayout_br1.addWidget(self.statToolbar_left)

        self.statFigure_right = Figure()
        self.statCanvas_right = FigureCanvas(self.statFigure_right)
        self.mainUI.verticalLayout_br2.addWidget(self.statCanvas_right)
        self.statToolbar_right = NavigationToolbar(self.statCanvas_right, self.mainwindow)
        self.mainUI.verticalLayout_br2.addWidget(self.statToolbar_right)

        self.addItems2StatsComboBox()

    def addItems2StatsComboBox(self):
        self.mainUI.comboBox_stats.insertItem(0, 'picked traces')
        self.mainUI.comboBox_stats.insertItem(1, 'mean SNR')
        self.mainUI.comboBox_stats.insertItem(2, 'median SNR')
        self.mainUI.comboBox_stats.insertItem(3, 'mean SPE')
        self.mainUI.comboBox_stats.insertItem(4, 'median SPE')
        self.enablePickedTools(False)

    def addItems2ShotsComboBox(self):
        shotnumbers = self.survey.data.keys()
        shotnumbers.sort()
        for index, shotnumber in enumerate(shotnumbers):
            self.mainUI.comboBox_shots.insertItem(index, 'Shot: %s'%shotnumber)
        self.mainUI.comboBox_shots.setMaxCount(len(shotnumbers))

    def increaseShotnumber(self):
        currentIndex = self.mainUI.comboBox_shots.currentIndex()
        maxindex = self.mainUI.comboBox_shots.maxCount() - 1
        if currentIndex == maxindex:
            self.mainUI.comboBox_shots.setCurrentIndex(0)
        else:
            self.mainUI.comboBox_shots.setCurrentIndex(currentIndex + 1)

    def decreaseShotnumber(self):
        currentIndex = self.mainUI.comboBox_shots.currentIndex()
        maxindex = self.mainUI.comboBox_shots.maxCount() - 1
        if currentIndex == 0:
            self.mainUI.comboBox_shots.setCurrentIndex(maxindex)
        else:
            self.mainUI.comboBox_shots.setCurrentIndex(currentIndex - 1)

        
    def plotShot(self):
        shotnumber = int(self.mainUI.comboBox_shots.currentText().split()[1])
        self.survey.data[shotnumber].matshow()

    def addArrayAxes(self):
        self.seisArrayAx = self.seisArrayFigure.add_subplot(111)

    def addSurfaceAxes(self):
        self.surfaceAx = self.surfaceFigure.add_subplot(111, projection = '3d')

    def addStatAxes(self):
        self.statAx_left = self.statFigure_left.add_subplot(111)
        self.statAx_right = self.statFigure_right.add_subplot(111)

    def enablePickedTools(self, bool, twoDim = False):
        self.mainUI.comboBox_stats.setEnabled(bool)
        self.statToolbar_left.setEnabled(bool)
        self.statToolbar_right.setEnabled(bool)
        if not twoDim:
            self.mainUI.comboBox_shots.setEnabled(bool)
            self.mainUI.shot_left.setEnabled(bool)
            self.mainUI.shot_right.setEnabled(bool)
            self.mainUI.plot_shot.setEnabled(bool)
        if bool == False:
            self.mainUI.comboBox_shots.clear()

    def replotArray(self):
        self.seisArrayFigure.clf()
        self.addArrayAxes()
        self.plotArray()
        self.seisArrayCanvas.draw()

    def replotSurface(self):
        self.surfaceFigure.clf()
        self.addSurfaceAxes()
        self.plotSurface()
        self.surfaceCanvas.draw()

    def plotArray(self):
        self.seisarray.plotArray2D(self.seisArrayAx, highlight_measured = True, plot_topo = True, twoDim = self.seisarray.twoDim)

    def plotSurface(self):
        if not self.seisarray.twoDim:
            self.seisarray.plotSurface3D(ax = self.surfaceAx, exag = True)
        self.seisarray.plotArray3D(ax = self.surfaceAx, legend = False, markersize = 3)

    def InitPickedWidgets(self):
        if self.checkPickState():
            surveyUtils.plotScatterStats4Receivers(self.survey, self.mainUI.comboBox_stats.currentText(),
                                                   self.statAx_left, twoDim = self.survey.twoDim)
            surveyUtils.plotScatterStats4Shots(self.survey, self.mainUI.comboBox_stats.currentText(),
                                               self.statAx_right, twoDim = self.survey.twoDim)
            self.addItems2ShotsComboBox()

    def refreshPickedWidgets(self):
        self.statFigure_left.clf()
        self.statFigure_right.clf()
        self.addStatAxes()
        self.InitPickedWidgets()
        self.statCanvas_left.draw()
        self.statCanvas_right.draw()

    def printSurveyTextbox(self, init = True):
        if init == True:
            surveytup = (0, 0, 0, 0)
        else:
            survey = self.survey
            nshots = len(survey.getShotlist())
            tt = survey.countAllTraces()
            pt = survey.countAllPickedTraces()
            rate = float(pt) / float(tt) * 100
            surveytup = (nshots, tt, pt, rate)
        surveyTitle = 'SURVEY:\n'
        surveyText = 'Number of Sources: %s\nTotal Traces: %s\nPicked Traces: %s (%4.2f%%)'%surveytup
        string = surveyTitle + surveyText
        self.mainUI.textBox_survey.setText(string)

    def printSeisArrayTextbox(self, init = True):
        if init == True:
            seistup = (0, 0, 0)
        else:
            seisarray = self.seisarray
            nshots = len(seisarray.getSourceCoordinates())
            nrec = len(seisarray.getReceiverCoordinates())
            nadd = len(seisarray.getMeasuredTopo())
            seistup = (nshots, nrec, nadd)
        seisArrayTitle = 'SEISARRAY:\n'
        seisArrayText = 'Sources: %s\nReceivers: %s\nAdditional Points:%s'%seistup
        string = seisArrayTitle + seisArrayText
        self.mainUI.textBox_seisarray.setText(string)

    def interpolate_receivers(self):
        if not self.checkSeisArrayState():
            self.printDialogMessage('No Seismic Array defined.')
            return
        self.seisarray.interpolateAll()
        self.refreshSeisArrayWidgets()

    def refreshSeisArrayWidgets(self):
        self.replotArray()
        self.replotSurface()
        self.printSeisArrayTextbox(init = False)

        
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
        self.printSurveyTextbox(init = False)
        print('Connected Seismic Array to active Survey object.')

    def startPicker(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        if self.checkPickState():
            if not self.continueDialogMessage('Survey already picked. Continue?'):
                return

        if self.autopicker is None:
            self.autopicker = Call_autopicker(self.mainwindow, self.survey)
        else:
            self.autopicker.start_dialog()
            self.autopicker.update_survey(self.survey)
        
        if self.autopicker.executed:
            self.setPickState(True)
            self.printSurveyTextbox(init = False)

    def startFMTOMO(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        if not self.checkPickState():
            self.printDialogMessage('Survey not picked.')
            return

        if self.fmtomo is None:
            self.fmtomo = Call_FMTOMO(self.mainwindow, self.survey)
        else:
            self.fmtomo.start_dialog()
            self.fmtomo.update_survey(self.survey)

        #if self.fmtomo.executed:
        

    def startVTKtools(self):
        if self.vtktools is None:
            self.vtktools = Call_VTK_dialog(self.mainwindow)
        else:
            self.vtktools.start_dialog()

    def postprocessing(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        self.postprocessing = Postprocessing(self.mainwindow, self.survey)
        #self.survey.plotAllPicks()
        #self.refreshPickedWidgets() # wait until finished

        
    def load_survey(self):
        if self.checkSurveyState():
            if not self.continueDialogExists('Survey'):
                return
        filename = openFile()
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

        filename = openFile()
        if filename is None:
            return
        try:
            seisarray = seismicArrayPreparation.SeisArray.from_pickle(filename)
        except:
            self.printDialogMessage('Could not load object %s.'%filename)
            return
        if not type(seisarray) == seismicArrayPreparation.SeisArray:
            self.printDialogMessage('Wrong input file of type %s, expected %s.'
                  %(type(seisarray), seismicArrayPreparation.SeisArray))
            return
        if disconnect:
            self.setConnected2SurveyState(False)
        self.seisarray = seisarray
        self.setSeisArrayState(True)

    def save_seisarray(self):
        if not self.checkSeisArrayState():
            self.printDialogMessage('No Seismic Array defined.')
            return
        filename = saveFile()
        if filename is None:
            return
        self.seisarray.saveSeisArray(filename)

    def save_survey(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        filename = saveFile()
        if filename is None:
            return
        self.survey.saveSurvey(filename)

    def setSurveyState(self, state):
        if state == True:
            self.mainUI.survey_active.setPixmap(self.applypixmap)
            self.printSurveyTextbox(init = False)
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
            self.refreshPickedWidgets()
            self.enablePickedTools(True, self.survey.twoDim)
            self.survey.picked = True
        elif state == True and self.checkSurveyState() is False:
            self.printDialogMessage('No Survey defined.')
            return
        elif state == False:
            self.mainUI.picked_active.setPixmap(self.cancelpixmap)
            if self.checkSurveyState():
                self.statFigure_left.clf()
                self.statFigure_right.clf()
                self.enablePickedTools(False)
                self.survey.picked = False

    def setSeisArrayState(self, state):
        if state == True:
            self.mainUI.seisarray_active.setPixmap(self.applypixmap)
            self.refreshSeisArrayWidgets()
            self.seisArrayToolbar.setEnabled(True)
        elif state == False:
            self.mainUI.seisarray_active.setPixmap(self.cancelpixmap)
            self.seisArrayToolbar.setEnabled(False)
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

    def exitApp(self):
        QtCore.QCoreApplication.instance().quit()


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.showMaximized()
    gui = gui_control()
    sys.exit(app.exec_())

