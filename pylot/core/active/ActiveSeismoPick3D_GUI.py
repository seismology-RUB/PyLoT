#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from PySide import QtCore, QtGui, QtCore
from asp3d_layout import *
from fmtomo_parameters_layout import *
from generate_survey_layout import *
from generate_seisarray_layout import *
from picking_parameters_layout import *
from pylot.core.active import activeSeismoPick, surveyUtils, fmtomoUtils, seismicArrayPreparation

class gui_control(object):
    def __init__(self):
        self.mainwindow = MainWindow
        self.mainUI = ui
        self.connectButtons()
        self.survey = None
        self.seisarray = None
        self.cancelpixmap = self.mainwindow.style().standardPixmap(QtGui.QStyle.SP_DialogCancelButton)
        self.applypixmap = self.mainwindow.style().standardPixmap(QtGui.QStyle.SP_DialogApplyButton)
        self.setInitStates()

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

    def gen_seisarray(self):
        if self.checkSeisArrayState():
            if not self.continueDialogExists('Seismic Array'):
                return
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
            self.setSeisArrayState(True)

    def gen_survey(self):
        if self.checkSurveyState():
            if not self.continueDialogExists('Survey'):
                return
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
            self.survey.setArtificialPick(0, 0) # artificial pick at source origin                                         
            surveyUtils.setDynamicFittedSNR(self.survey.getShotDict())
            self.setSurveyState(True)

    def interpolate_receivers(self):
        if not self.seisarray:
            print('No Seismic Array defined.')
        self.seisarray.interpolateAll()

    def getPickParameters(self, ui, Picking_parameters):
        if Picking_parameters.exec_():
            ncores = int(ui.lineEdit_ncores.text())
            vmin = float(ui.lineEdit_vmin.text())
            vmax = float(ui.lineEdit_vmax.text())
            folm = float(ui.lineEdit_folm.text())
            AIC = ui.checkBox.isChecked()
            aicwindow = [int(val) for val in ui.lineEdit_aicwindow.text().split(',')]
            return ncores, vmin, vmax, folm, AIC, tuple(aicwindow)

    def connect2Survey(self):
        if not self.checkSurveyState():
            print('Survey not defined.')
            return
        if not self.checkSeisArrayState():
            print('Got no Seismic Array.')
            return
        if self.checkConnected2SurveyState():
            if self.continueDialogMessage('Existing Survey already got Seismic Array object. Continue?'):
                pass
        self.survey.seisarray = self.seisarray
        self.setConnected2SurveyState(True)
        print('Connected Seismic Array to active Survey object.')

    def callPicker(self):
        if not self.checkSurveyState():
            print('Survey not defined.')
            return
        if self.checkPickState():
            if not self.continueDialogMessage('Survey already picked. Continue?'):
                return
        Picking_parameters = QtGui.QDialog(self.mainwindow)
        ui = Ui_Picking_parameters()
        ui.setupUi(Picking_parameters)
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
            print('Survey not defined.')
            return
        fmtomo_parameters = QtGui.QDialog(self.mainwindow)
        ui = Ui_fmtomo_parameters()
        ui.setupUi(fmtomo_parameters)                
        self.fmtomo_parameters_ui = ui
        self.connectButtons_startFMTOMO()
        self.getFMTOMOparameters(ui, fmtomo_parameters)

    def getFMTOMOparameters(self, ui, fmtomo_parameters):
        if fmtomo_parameters.exec_():
            fmtomo_dir = ui.fmtomo_dir.text()
            nIter = int(ui.nIter.value())
            nproc = int(ui.nproc.text())
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
        QtCore.QObject.connect(self.fmtomo_parameters_ui.browse_seisarray, QtCore.SIGNAL("clicked()"), self.chooseSeisarray)

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
            print('Survey not defined.')
            return
        self.survey.plotAllPicks()
        
    def load_survey(self):
        if self.checkSurveyState():
            if not self.continueDialogExists('Survey'):
                return
        filename = self.openFile()
        if filename is None:
            return
        self.survey = activeSeismoPick.Survey.from_pickle(filename)
        self.setSurveyState(True)
        if self.survey.picked:
            self.setPickState(True)
        else:
            self.setPickState(False)
        if self.survey.seisarray != None:
            self.seisarray = self.survey.seisarray
            self.setConnected2SurveyState(True)
            self.setSeisArrayState(True)
            print('Loaded Survey with active Seismic Array.')

    def load_seisarray(self):
        if self.checkSeisArrayState():
            if not self.continueDialogExists('Seismic Array'):
                return
        filename = self.openFile()
        if filename is None:
            return
        self.seisarray = seismicArrayPreparation.SeisArray.from_pickle(filename)
        self.setSeisArrayState(True)

    def save_seisarray(self):
        if not self.checkSeisArrayState():
            print('No Seismic Array defined.')
            return
        filename = self.saveFile()
        if filename is None:
            return
        self.seisarray.saveSeisArray(filename)

    def save_survey(self):
        if not self.checkSurveyState():
            print('No Survey defined.')
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
            self.survey.picked = True
        elif state == True and self.checkSurveyState() is False:
            print('No Survey defined.')
            return
        elif state == False:
            self.mainUI.picked_active.setPixmap(self.cancelpixmap)
            if self.checkSurveyState():
                self.survey.picked = False

    def setSeisArrayState(self, state):
        if state == True:
            self.mainUI.seisarray_active.setPixmap(self.applypixmap)
        elif state == False:
            self.mainUI.seisarray_active.setPixmap(self.cancelpixmap)

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
            print ('No Survey defined.')
            return
        return self.survey.picked

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

    def connectButtons_gen_seisarray(self):
        QtCore.QObject.connect(self.gen_new_seisarray.pushButton_rec, QtCore.SIGNAL("clicked()"), self.chooseMeasuredRec)
        QtCore.QObject.connect(self.gen_new_seisarray.pushButton_src, QtCore.SIGNAL("clicked()"), self.chooseMeasuredSrc)
        QtCore.QObject.connect(self.gen_new_seisarray.pushButton_obs, QtCore.SIGNAL("clicked()"), self.chooseMeasuredPts)

    def chooseMeasuredSrc(self):
        self.gen_new_seisarray.lineEdit_src.setText(self.openFile())

    def chooseMeasuredRec(self):
        self.gen_new_seisarray.lineEdit_rec.setText(self.openFile())

    def chooseMeasuredPts(self):
        self.gen_new_seisarray.lineEdit_pts.setText(self.browseDir())

    def chooseSourcefile(self):
        self.gen_new_survey.lineEdit_src.setText(self.openFile())

    def chooseReceiverfile(self):
        self.gen_new_survey.lineEdit_rec.setText(self.openFile())

    def chooseObsdir(self):
        self.gen_new_survey.lineEdit_obs.setText(self.browseDir())

    def openFile(self):
        dialog = QtGui.QFileDialog()
        filename = dialog.getOpenFileName()
        if len(filename[0]) > 0:
            return filename[0]

    def saveFile(self):
        dialog = QtGui.QFileDialog()
        filename = dialog.getSaveFileName()
        if len(filename[0]) > 0:
            return filename[0]

    def browseDir(self):
        dialog = QtGui.QFileDialog()
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


