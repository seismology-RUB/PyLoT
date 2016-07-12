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
from vtk_tools_layout import *
from pylot.core.active import activeSeismoPick, surveyUtils, fmtomoUtils, seismicArrayPreparation

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
        QtCore.QObject.connect(self.mainUI.vtk_tools, QtCore.SIGNAL("clicked()"), self.startVTKtools)
        QtCore.QObject.connect(self.mainUI.comboBox_stats, QtCore.SIGNAL("activated(int)"), self.refreshPickedWidgets)
        QtCore.QObject.connect(self.mainUI.shot_left, QtCore.SIGNAL("clicked()"), self.decreaseShotnumber)
        QtCore.QObject.connect(self.mainUI.shot_right, QtCore.SIGNAL("clicked()"), self.increaseShotnumber)
        QtCore.QObject.connect(self.mainUI.plot_shot, QtCore.SIGNAL("clicked()"), self.plotShot)

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
            if self.gen_new_seisarray.radioButton_interpolatable.isChecked():
                self.seisarray = seismicArrayPreparation.SeisArray(recfile, True)
            elif self.gen_new_seisarray.radioButton_normal.isChecked():
                self.seisarray = seismicArrayPreparation.SeisArray(recfile, False)
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
            fstart = self.gen_new_survey_min.fstart.text()
            fend = self.gen_new_survey_min.fend.text()
            self.survey = activeSeismoPick.Survey(obsdir, seisArray = self.seisarray,
                                                  useDefaultParas = True, fstart = fstart,
                                                  fend = fend)
            self.setConnected2SurveyState(True)
            self.setPickState(False)
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
            fstart = self.gen_new_survey.fstart.text()
            fend = self.gen_new_survey.fend.text()
            self.survey = activeSeismoPick.Survey(obsdir, srcfile, recfile,
                                                  useDefaultParas = True,
                                                  fstart = fstart, fend = fend)
            self.setConnected2SurveyState(False)
            return True

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

    def refreshPickedWidgets(self):
        self.statFigure_left.clf()
        self.statFigure_right.clf()
        self.addStatAxes()
        self.InitPickedWidgets()
        self.statCanvas_left.draw()
        self.statCanvas_right.draw()

    def InitPickedWidgets(self):
        if self.checkPickState():
            surveyUtils.plotScatterStats4Receivers(self.survey, self.mainUI.comboBox_stats.currentText(),
                                                   self.statAx_left, twoDim = self.survey.twoDim)
            surveyUtils.plotScatterStats4Shots(self.survey, self.mainUI.comboBox_stats.currentText(),
                                               self.statAx_right, twoDim = self.survey.twoDim)
            self.addItems2ShotsComboBox()

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
        
    def getPickParameters(self, ui, Picking_parameters):
        if Picking_parameters.exec_():
            ncores = int(ui.ncores.value())
            vmin = float(ui.lineEdit_vmin.text())
            vmax = float(ui.lineEdit_vmax.text())
            folm = float(ui.slider_folm.value())/100.
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
        self.printSurveyTextbox(init = False)
        print('Connected Seismic Array to active Survey object.')

    def getMaxCPU(self):
        import multiprocessing
        return multiprocessing.cpu_count()

    def refreshFolm(self):
        self.picker_ui.label_folm.setText('%s %%'%self.picker_ui.slider_folm.value())

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
        self.picker_ui = ui
        QtCore.QObject.connect(self.picker_ui.slider_folm, QtCore.SIGNAL("valueChanged(int)"), self.refreshFolm)
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
        self.printSurveyTextbox(init = False)

    def startFMTOMO(self):
        if not self.checkSurveyState():
            self.printDialogMessage('No Survey defined.')
            return
        if not self.checkPickState():
            self.printDialogMessage('Survey not picked.')
            return
        fmtomo_parameters = QtGui.QDialog(self.mainwindow)
        ui = Ui_fmtomo_parameters()
        ui.setupUi(fmtomo_parameters)                
        ui.nproc.setMaximum(self.getMaxCPU())

        self.fmtomo_parameters_ui = ui
        self.connectButtons_startFMTOMO()
        self.getFMTOMOparameters(ui, fmtomo_parameters)

    def startVTKtools(self):
        vtk_tools = QtGui.QDialog(self.mainwindow)
        ui = Ui_vtk_tools()
        ui.setupUi(vtk_tools)

        self.vtk_tools_ui = ui
        self.connectButtons_vtk_tools()
        self.openVTKdialog(ui, vtk_tools)

    def openVTKdialog(self, ui, vtk_tools):
        vtk_tools.exec_()

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
            if self.seisarray.twoDim:
                interpolationMethod = 'nearest'
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

    def connectButtons_vtk_tools(self):
        QtCore.QObject.connect(self.vtk_tools_ui.pushButton_vg, QtCore.SIGNAL("clicked()"), self.chooseVgrid)
        QtCore.QObject.connect(self.vtk_tools_ui.pushButton_vgref, QtCore.SIGNAL("clicked()"), self.chooseVgridref)
        QtCore.QObject.connect(self.vtk_tools_ui.pushButton_rays, QtCore.SIGNAL("clicked()"), self.chooseRaysIn)
        QtCore.QObject.connect(self.vtk_tools_ui.pushButton_raysout, QtCore.SIGNAL("clicked()"), self.chooseRaysOutDir)
        QtCore.QObject.connect(self.vtk_tools_ui.pushButton_vtkout, QtCore.SIGNAL("clicked()"), self.newFileVTK)
        QtCore.QObject.connect(self.vtk_tools_ui.pushButton_parav, QtCore.SIGNAL("clicked()"), self.openFileParaview)
        QtCore.QObject.connect(self.vtk_tools_ui.start_vg, QtCore.SIGNAL("clicked()"), self.startvgvtk)
        QtCore.QObject.connect(self.vtk_tools_ui.start_rays, QtCore.SIGNAL("clicked()"), self.startraysvtk)
        QtCore.QObject.connect(self.vtk_tools_ui.radioButton_rel, QtCore.SIGNAL("clicked()"), self.activateVgref)
        QtCore.QObject.connect(self.vtk_tools_ui.radioButton_abs, QtCore.SIGNAL("clicked()"), self.deactivateVgref)

    def openFileParaview(self):
        os.system('paraview %s &'%self.vtk_tools_ui.lineEdit_vgout.text())

    def activateVgref(self):
        self.vtk_tools_ui.lineEdit_vgref.setEnabled(True)
        self.vtk_tools_ui.pushButton_vgref.setEnabled(True)

    def deactivateVgref(self):
        self.vtk_tools_ui.lineEdit_vgref.setEnabled(False)
        self.vtk_tools_ui.pushButton_vgref.setEnabled(False)

    def checkVgStartButton(self):
        ui = self.vtk_tools_ui
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
        ui = self.vtk_tools_ui
        if ui.lineEdit_rays.text() != '' and ui.lineEdit_raysout.text() != '':
            ui.start_rays.setEnabled(True)
        else:
            ui.start_rays.setEnabled(False)

    def chooseVgrid(self):
        self.vtk_tools_ui.lineEdit_vg.setText(self.openFile())
        self.checkVgStartButton()

    def chooseVgridref(self):
        self.vtk_tools_ui.lineEdit_vgref.setText(self.openFile())
        self.checkVgStartButton()

    def chooseRaysIn(self):
        self.vtk_tools_ui.lineEdit_rays.setText(self.openFile())
        self.checkRaysStartButton()

    def chooseRaysOutDir(self):
        self.vtk_tools_ui.lineEdit_raysout.setText(self.browseDir())
        self.checkRaysStartButton()

    def startvgvtk(self):
        ui = self.vtk_tools_ui
        if ui.lineEdit_vgout.text() == '':
            if not self.printDialogMessage('Please specify output filename.'):
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
        ui = self.vtk_tools_ui
        fmtomoUtils.rays2VTK(ui.lineEdit_rays.text(), ui.lineEdit_raysout.text())

    def newFileVTK(self):
        self.vtk_tools_ui.lineEdit_vgout.setText(self.saveFile())

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
        self.refreshPickedWidgets() # wait until finished
        
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


