# -*- coding: utf-8 -*-

import sys
from PySide import QtCore, QtGui
from pylot.core.active import activeSeismoPick

class Ui_ActiveSeismoPick3D(object):
    def __init__(self):
        self.survey = None
        
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(550, 350)
        
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.addBrowseButtons()
        self.addButtons()
        self.addLineEdits()
        self.addLabels()
        
        MainWindow.setCentralWidget(self.centralwidget)

        self.setMenubar(MainWindow)
        
        self.menuPreferences = QtGui.QMenu(self.menubar)
        self.menuPreferences.setObjectName("menuPreferences")
        MainWindow.setMenuBar(self.menubar)

        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.menubar.addAction(self.menuPreferences.menuAction())

        self.retranslateUi(MainWindow)
        self.connectButtons()

    def setMenubar(self, window):
        self.menubar = QtGui.QMenuBar(window)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 603, 21))
        self.menubar.setObjectName("menubar")
        
    def connectButtons(self):
        QtCore.QObject.connect(self.browseButton_rec, QtCore.SIGNAL("clicked()"), self.chooseReceiverfile)
        QtCore.QObject.connect(self.browseButton_src, QtCore.SIGNAL("clicked()"), self.chooseSourcefile)
        QtCore.QObject.connect(self.browseButton_obsdir, QtCore.SIGNAL("clicked()"), self.chooseObsdir)
        QtCore.QObject.connect(self.start_picking, QtCore.SIGNAL("clicked()"), self.callPicker)
        QtCore.QObject.connect(self.gen_survey, QtCore.SIGNAL("clicked()"), self.generateSurvey)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        
    def addLabels(self):
        self.label = QtGui.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(40, 70, 131, 20))
        self.label.setObjectName("label")
        
        self.label_2 = QtGui.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(40, 110, 131, 20))
        self.label_2.setObjectName("label_2")
        
        self.label_3 = QtGui.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(40, 150, 131, 20))
        self.label_3.setObjectName("label_3")

    def addLineEdits(self):
        self.receiverfile_lineEdit = QtGui.QLineEdit(self.centralwidget)
        self.receiverfile_lineEdit.setGeometry(QtCore.QRect(192, 70, 231, 20))
        self.receiverfile_lineEdit.setObjectName("receiverfile_lineEdit")

        self.sourcefile_lineEdit = QtGui.QLineEdit(self.centralwidget)
        self.sourcefile_lineEdit.setGeometry(QtCore.QRect(192, 110, 231, 20))
        self.sourcefile_lineEdit.setObjectName("sourcefile_lineEdit")

        self.obsdir_lineEdit = QtGui.QLineEdit(self.centralwidget)
        self.obsdir_lineEdit.setGeometry(QtCore.QRect(192, 150, 231, 20))
        self.obsdir_lineEdit.setText("")
        self.obsdir_lineEdit.setObjectName("obsdir_lineEdit")

    def addBrowseButtons(self):
        self.browseButton_rec = QtGui.QPushButton(self.centralwidget)
        self.browseButton_rec.setGeometry(QtCore.QRect(440, 70, 75, 23))
        self.browseButton_rec.setObjectName("browseButton_rec")
        
        self.browseButton_src = QtGui.QPushButton(self.centralwidget)
        self.browseButton_src.setGeometry(QtCore.QRect(440, 110, 75, 23))
        self.browseButton_src.setObjectName("browseButton_src")

        self.browseButton_obsdir = QtGui.QPushButton(self.centralwidget)
        self.browseButton_obsdir.setGeometry(QtCore.QRect(440, 150, 75, 23))
        self.browseButton_obsdir.setObjectName("browseButton_obsdir")

    def addButtons(self):
        self.gen_survey = QtGui.QPushButton(self.centralwidget)
        self.gen_survey.setGeometry(QtCore.QRect(80, 230, 61, 61))
        self.gen_survey.setObjectName("gen_survey")

        self.start_picking = QtGui.QPushButton(self.centralwidget)
        self.start_picking.setGeometry(QtCore.QRect(160, 230, 61, 61))
        self.start_picking.setObjectName("start_picking")

    def addProgressBar(self):
        self.progressBar = QtGui.QProgressBar(self.centralwidget)
        self.progressBar.setGeometry(QtCore.QRect(470, 280, 118, 23))
        self.progressBar.setObjectName("progressBar")

    def updateProgressBar(self):
        self.progressBar.setProperty("value", 24)

    def browseFile(self):
        dialog = QtGui.QFileDialog()
        filename = dialog.getOpenFileName()
        return filename[0]

    def browseDir(self):
        dialog = QtGui.QFileDialog()
        directory = dialog.getExistingDirectory()
        return directory

    def chooseSourcefile(self):
        self.sourcefile_lineEdit.setText(self.browseFile())

    def chooseReceiverfile(self):
        self.receiverfile_lineEdit.setText(self.browseFile())

    def chooseObsdir(self):
        self.obsdir_lineEdit.setText(self.browseDir())
    
    def getSourcefile(self):
        return self.sourcefile_lineEdit.text()

    def getReceiverfile(self):
        return self.receiverfile_lineEdit.text()

    def getObsdir(self):
        return self.obsdir_lineEdit.text()
    
    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "ActiveSeismoPick3D", None, QtGui.QApplication.UnicodeUTF8))
        self.browseButton_rec.setText(QtGui.QApplication.translate("MainWindow", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.browseButton_src.setText(QtGui.QApplication.translate("MainWindow", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Receiver File", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Source File", None, QtGui.QApplication.UnicodeUTF8))
        self.start_picking.setText(QtGui.QApplication.translate("MainWindow", "Start\n"
"Picking", None, QtGui.QApplication.UnicodeUTF8))
        self.browseButton_obsdir.setText(QtGui.QApplication.translate("MainWindow", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.gen_survey.setText(QtGui.QApplication.translate("MainWindow", "Gemerate\nSurvey", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Seismogram Directory", None, QtGui.QApplication.UnicodeUTF8))
        self.menuPreferences.setTitle(QtGui.QApplication.translate("MainWindow", "Preferences", None, QtGui.QApplication.UnicodeUTF8))

    def generateSurvey(self):
        obsdir = self.getObsdir()
        self.survey = activeSeismoPick.Survey(self.getObsdir(), self.getSourcefile(), self.getReceiverfile(),
                                         useDefaultParas = True)
                                         
        
    def callPicker(self):
        Picking_parameters = QtGui.QDialog(self.centralwidget)
        ui = Ui_Picking_parameters()
        ui.setupUi(Picking_parameters)
        ncores, vmin, vmax, folm, AIC, aicwindow = ui.getParameters(Picking_parameters)
        if AIC == True:
            HosAic = 'aic'
        else:
            HosAic = 'hos'
        if self.survey is None:
            print('Survey not defined.')
            return
                
        self.survey.pickAllShots(vmin = vmin, vmax = vmax,
                                 folm = folm, HosAic = HosAic,
                                 aicwindow = aicwindow, cores = ncores)
        
    
class Ui_Picking_parameters(object):
    def setupUi(self, Picking_parameters):
        Picking_parameters.setObjectName("Picking_parameters")
        Picking_parameters.resize(300, 300)
        self.gridLayoutWidget = QtGui.QWidget(Picking_parameters)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(20, 20, 250, 211))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtGui.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.label_cores = QtGui.QLabel(self.gridLayoutWidget)
        self.label_cores.setObjectName("label_cores")
        self.gridLayout.addWidget(self.label_cores, 0, 0, 1, 1)
        self.label_vmax = QtGui.QLabel(self.gridLayoutWidget)
        self.label_vmax.setObjectName("label_vmax")
        self.gridLayout.addWidget(self.label_vmax, 2, 0, 1, 1)
        self.label_vmin = QtGui.QLabel(self.gridLayoutWidget)
        self.label_vmin.setObjectName("label_vmin")
        self.gridLayout.addWidget(self.label_vmin, 1, 0, 1, 1)
        self.lineEdit_ncores = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_ncores.setObjectName("lineEdit_ncores")
        self.gridLayout.addWidget(self.lineEdit_ncores, 0, 1, 1, 1)
        self.lineEdit_vmin = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_vmin.setObjectName("lineEdit_vmin")
        self.gridLayout.addWidget(self.lineEdit_vmin, 1, 1, 1, 1)
        self.lineEdit_vmax = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_vmax.setObjectName("lineEdit_vmax")
        self.gridLayout.addWidget(self.lineEdit_vmax, 2, 1, 1, 1)
        self.checkBox = QtGui.QCheckBox(self.gridLayoutWidget)
        self.checkBox.setObjectName("checkBox")
        self.gridLayout.addWidget(self.checkBox, 4, 1, 1, 1)
        self.label_folm = QtGui.QLabel(self.gridLayoutWidget)
        self.label_folm.setObjectName("label_folm")
        self.gridLayout.addWidget(self.label_folm, 3, 0, 1, 1)
        self.label_aic = QtGui.QLabel(self.gridLayoutWidget)
        self.label_aic.setObjectName("label_aic")
        self.gridLayout.addWidget(self.label_aic, 4, 0, 1, 1)
        self.lineEdit_folm = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_folm.setObjectName("lineEdit_folm")
        self.gridLayout.addWidget(self.lineEdit_folm, 3, 1, 1, 1)
        self.label_aicwindow = QtGui.QLabel(self.gridLayoutWidget)
        self.label_aicwindow.setObjectName("label_aicwindow")
        self.gridLayout.addWidget(self.label_aicwindow, 5, 0, 1, 1)
        self.lineEdit_aicwindow = QtGui.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_aicwindow.setObjectName("lineEdit_aicwindow")
        self.gridLayout.addWidget(self.lineEdit_aicwindow, 5, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(Picking_parameters)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.buttonBox.setGeometry(QtCore.QRect(10, 240, 250, 32))

        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), Picking_parameters.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), Picking_parameters.reject)
        self.retranslateUi(Picking_parameters)
        QtCore.QMetaObject.connectSlotsByName(Picking_parameters)

    def getParameters(self, Picking_parameters):
        if Picking_parameters.exec_():
            ncores = int(self.lineEdit_ncores.text())
            vmin = float(self.lineEdit_vmin.text())
            vmax = float(self.lineEdit_vmax.text())
            folm = float(self.lineEdit_folm.text())
            AIC = self.checkBox.isChecked()
            aicwindow = [float(val) for val in self.lineEdit_aicwindow.text().split(',')]
            
            return ncores, vmin, vmax, folm, AIC, tuple(aicwindow)
            
    def retranslateUi(self, Picking_parameters):
        Picking_parameters.setWindowTitle(QtGui.QApplication.translate("Picking_parameters", "Picking_parameters", None, QtGui.QApplication.UnicodeUTF8))
        self.label_cores.setText(QtGui.QApplication.translate("Picking_parameters", "Number of cores", None, QtGui.QApplication.UnicodeUTF8))
        self.label_vmax.setText(QtGui.QApplication.translate("Picking_parameters", "Vmax (default = 5000m/s)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_vmin.setText(QtGui.QApplication.translate("Picking_parameters", "Vmin (default = 333 m/s)", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_ncores.setText(QtGui.QApplication.translate("Picking_parameters", "1", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_vmin.setText(QtGui.QApplication.translate("Picking_parameters", "333", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_vmax.setText(QtGui.QApplication.translate("Picking_parameters", "5000", None, QtGui.QApplication.UnicodeUTF8))
        self.label_folm.setText(QtGui.QApplication.translate("Picking_parameters", "Fraction of local maximum\n"
"(default = 0.6)", None, QtGui.QApplication.UnicodeUTF8))
        self.label_aic.setText(QtGui.QApplication.translate("Picking_parameters", "AIC", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_folm.setText(QtGui.QApplication.translate("Picking_parameters", "0.6", None, QtGui.QApplication.UnicodeUTF8))
        self.label_aicwindow.setText(QtGui.QApplication.translate("Picking_parameters", "AIC window (only if AIC\n"
" is checked)", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_aicwindow.setText(QtGui.QApplication.translate("Picking_parameters", "15, 0", None, QtGui.QApplication.UnicodeUTF8))

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_ActiveSeismoPick3D()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

