#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
from PySide import QtCore, QtGui, QtCore
from asp3d_layout import *
from pylot.core.active import activeSeismoPick, surveyUtils, fmtomoUtils

class gui_control(object):
    def __init__(self):
        self.mainwindow = MainWindow
        self.mainUI = ui
        self.connectButtons()
        self.survey = None
        
    def connectButtons(self):
        QtCore.QObject.connect(self.mainUI.gen_new_survey, QtCore.SIGNAL("clicked()"), self.gen_survey)
        QtCore.QObject.connect(self.mainUI.load_survey, QtCore.SIGNAL("clicked()"), self.load_survey)
        QtCore.QObject.connect(self.mainUI.save_survey, QtCore.SIGNAL("clicked()"), self.save_survey)
        QtCore.QObject.connect(self.mainUI.picker, QtCore.SIGNAL("clicked()"), self.callPicker)
        QtCore.QObject.connect(self.mainUI.postprocessing, QtCore.SIGNAL("clicked()"), self.postprocessing)
        QtCore.QObject.connect(self.mainUI.fmtomo, QtCore.SIGNAL("clicked()"), self.startFMTOMO)

    def gen_survey(self):
        if self.checkSurveyState():
            if not self.continueSurveyDialog():
                return
        qdialog = QtGui.QDialog(self.mainwindow)
        ui = Ui_Dialog()
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
            self.setSurveyState('active')

    def getPickParameters(self, ui, Picking_parameters):
        if Picking_parameters.exec_():
            ncores = int(ui.lineEdit_ncores.text())
            vmin = float(ui.lineEdit_vmin.text())
            vmax = float(ui.lineEdit_vmax.text())
            folm = float(ui.lineEdit_folm.text())
            AIC = ui.checkBox.isChecked()
            aicwindow = [int(val) for val in ui.lineEdit_aicwindow.text().split(',')]
            return ncores, vmin, vmax, folm, AIC, tuple(aicwindow)
            
    def callPicker(self):
        if self.survey is None:
            print('Survey not defined.')
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
        self.setPickState('active')

    def startFMTOMO(self):
        if self.survey is None:
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
            seisarray_loc = ui.seisarray.text()

            if not os.path.isdir(picks_dir):
                err = os.mkdir(picks_dir)

            self.survey.exportFMTOMO(picks_dir)
            self.survey.loadArrayFromPickle(seisarray_loc)

            cwd = os.getcwd()
            interpolationMethod = 'linear'
            os.chdir(simuldir)
            self.survey.seisArray.generateFMTOMOinputFromArray(propgrid, vgrid, (bbot, btop), cushionfactor,
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
        self.fmtomo_parameters_ui.customgrid.setText(self.browseFile())

    def chooseSimuldir(self):
        self.fmtomo_parameters_ui.simuldir.setText(self.browseDir())

    def chooseSeisarray(self):
        self.fmtomo_parameters_ui.seisarray.setText(self.browseFile())

    def postprocessing(self):
        if self.survey is None:
            print('Survey not defined.')
            return
        self.survey.plotAllPicks()
        
    def load_survey(self):
        if self.checkSurveyState():
            if not self.continueSurveyDialog():
                return
        filename = self.browseFile()
        if filename is None:
            return
        self.survey = activeSeismoPick.Survey.from_pickle(filename)
        self.setSurveyState('active')
        if self.survey.picked:
            self.setPickState('active')

    def save_survey(self):
        if not self.checkSurveyState():
            print('No Survey defined.')
            return
        if self.checkSurveyState:
            filename = self.browseFile()
            if filename is None:
                return
            self.survey.saveSurvey(filename)
        else:
            print('No active Survey.')

    def setSurveyState(self, state):
        if state == 'active':
            self.mainUI.survey_active.setCheckable(True)
            self.mainUI.survey_active.setChecked(True)
            self.mainUI.survey_active.setCheckable(False)
        elif state == 'inactive':
            self.mainUI.survey_active.setCheckable(True)
            self.mainUI.survey_active.setChecked(False)
            self.mainUI.survey_active.setCheckable(False)

    def checkSurveyState(self):
        if self.mainUI.survey_active.checkState() == QtCore.Qt.CheckState.Checked:
            return True
        if self.mainUI.survey_active.checkState() == QtCore.Qt.CheckState.Unchecked:
            return False

    def setPickState(self, state):
        if state == 'active' and self.checkSurveyState():
            self.mainUI.picked_active.setCheckable(True)
            self.mainUI.picked_active.setChecked(True)
            self.mainUI.picked_active.setCheckable(True)
        elif self.checkSurveyState() is False:
            print('No Survey defined.')
            return
        elif state == 'inactive':
            self.mainUI.picked_active.setCheckable(True)
            self.mainUI.picked_active.setChecked(False)
            self.mainUI.picked_active.setCheckable(True)

    def checkPickState(self):
        if self.mainUI.picked_active.checkState() == QtCore.Qt.CheckState.Checked:
            return True
        if self.mainUI.picked_active.checkState() == QtCore.Qt.CheckState.Unchecked:
            return False

    def continueSurveyDialog(self):
        qmb = QtGui.QMessageBox()
        qmb.setText('Survey already exists. Overwrite?')
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

    def chooseSourcefile(self):
        self.gen_new_survey.lineEdit_src.setText(self.browseFile())

    def chooseReceiverfile(self):
        self.gen_new_survey.lineEdit_rec.setText(self.browseFile())

    def chooseObsdir(self):
        self.gen_new_survey.lineEdit_obs.setText(self.browseDir())

    def browseFile(self):
        dialog = QtGui.QFileDialog()
        filename = dialog.getOpenFileName()
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


