# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'generate_survey_layout_minimal.ui'
#
# Created: Wed Jul  6 11:39:58 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_generate_survey_minimal(object):
    def setupUi(self, generate_survey_minimal):
        generate_survey_minimal.setObjectName("generate_survey_minimal")
        generate_survey_minimal.resize(325, 83)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../asp3d_icon.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        generate_survey_minimal.setWindowIcon(icon)
        self.verticalLayout = QtGui.QVBoxLayout(generate_survey_minimal)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.lineEdit_obs = QtGui.QLineEdit(generate_survey_minimal)
        self.lineEdit_obs.setObjectName("lineEdit_obs")
        self.gridLayout.addWidget(self.lineEdit_obs, 0, 1, 1, 1)
        self.label_obs = QtGui.QLabel(generate_survey_minimal)
        self.label_obs.setObjectName("label_obs")
        self.gridLayout.addWidget(self.label_obs, 0, 0, 1, 1)
        self.pushButton_obs = QtGui.QPushButton(generate_survey_minimal)
        self.pushButton_obs.setObjectName("pushButton_obs")
        self.gridLayout.addWidget(self.pushButton_obs, 0, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.buttonBox = QtGui.QDialogButtonBox(generate_survey_minimal)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(generate_survey_minimal)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), generate_survey_minimal.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), generate_survey_minimal.reject)
        QtCore.QMetaObject.connectSlotsByName(generate_survey_minimal)

    def retranslateUi(self, generate_survey_minimal):
        generate_survey_minimal.setWindowTitle(QtGui.QApplication.translate("generate_survey_minimal", "Generate new Survey", None, QtGui.QApplication.UnicodeUTF8))
        self.label_obs.setToolTip(QtGui.QApplication.translate("generate_survey_minimal", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Specifiy directory containing seismograms for each shot.</p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Currently in the format SEGY with each file named \'shotnumber*_pickle.dat\'.</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Shot number 100 containing seismograms for all traces with the name:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">100_pickle.dat</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_obs.setText(QtGui.QApplication.translate("generate_survey_minimal", "Seismogram\n"
"Directory [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_obs.setText(QtGui.QApplication.translate("generate_survey_minimal", "Browse", None, QtGui.QApplication.UnicodeUTF8))

