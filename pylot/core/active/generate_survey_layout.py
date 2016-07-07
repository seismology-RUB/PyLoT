# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'generate_survey_layout.ui'
#
# Created: Thu Jul  7 14:25:26 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_generate_survey(object):
    def setupUi(self, generate_survey):
        generate_survey.setObjectName("generate_survey")
        generate_survey.resize(380, 160)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../asp3d_icon.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        generate_survey.setWindowIcon(icon)
        self.verticalLayout = QtGui.QVBoxLayout(generate_survey)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.lineEdit_rec = QtGui.QLineEdit(generate_survey)
        self.lineEdit_rec.setObjectName("lineEdit_rec")
        self.gridLayout.addWidget(self.lineEdit_rec, 0, 1, 1, 1)
        self.pushButton_rec = QtGui.QPushButton(generate_survey)
        self.pushButton_rec.setObjectName("pushButton_rec")
        self.gridLayout.addWidget(self.pushButton_rec, 0, 2, 1, 1)
        self.label_rec = QtGui.QLabel(generate_survey)
        self.label_rec.setObjectName("label_rec")
        self.gridLayout.addWidget(self.label_rec, 0, 0, 1, 1)
        self.lineEdit_obs = QtGui.QLineEdit(generate_survey)
        self.lineEdit_obs.setObjectName("lineEdit_obs")
        self.gridLayout.addWidget(self.lineEdit_obs, 2, 1, 1, 1)
        self.label_obs = QtGui.QLabel(generate_survey)
        self.label_obs.setObjectName("label_obs")
        self.gridLayout.addWidget(self.label_obs, 2, 0, 1, 1)
        self.pushButton_obs = QtGui.QPushButton(generate_survey)
        self.pushButton_obs.setObjectName("pushButton_obs")
        self.gridLayout.addWidget(self.pushButton_obs, 2, 2, 1, 1)
        self.label_src = QtGui.QLabel(generate_survey)
        self.label_src.setObjectName("label_src")
        self.gridLayout.addWidget(self.label_src, 1, 0, 1, 1)
        self.lineEdit_src = QtGui.QLineEdit(generate_survey)
        self.lineEdit_src.setObjectName("lineEdit_src")
        self.gridLayout.addWidget(self.lineEdit_src, 1, 1, 1, 1)
        self.pushButton_src = QtGui.QPushButton(generate_survey)
        self.pushButton_src.setObjectName("pushButton_src")
        self.gridLayout.addWidget(self.pushButton_src, 1, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.buttonBox = QtGui.QDialogButtonBox(generate_survey)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(generate_survey)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), generate_survey.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), generate_survey.reject)
        QtCore.QMetaObject.connectSlotsByName(generate_survey)

    def retranslateUi(self, generate_survey):
        generate_survey.setWindowTitle(QtGui.QApplication.translate("generate_survey", "Generate new Survey", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_rec.setText(QtGui.QApplication.translate("generate_survey", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_rec.setToolTip(QtGui.QApplication.translate("generate_survey", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Load receiver input file. The input file must be in the following format:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Containing in each line, separated by spaces:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">[trace ID (int)] [X (float)]  [Y (float)]  [Z (float)]</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Trace ID 100 with the coordinates (12.3 [m], 100.5 [m], 20.3 [m]).</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">100  12.3  100.5  20.3</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_rec.setText(QtGui.QApplication.translate("generate_survey", "Receiver\n"
"File [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.label_obs.setToolTip(QtGui.QApplication.translate("generate_survey", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
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
        self.label_obs.setText(QtGui.QApplication.translate("generate_survey", "Seismogram\n"
"Directory [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_obs.setText(QtGui.QApplication.translate("generate_survey", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setToolTip(QtGui.QApplication.translate("generate_survey", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Load sources input file. The input file must be in the following format:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Containing in each line, separated by spaces:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">[trace ID (int)] [X (float)]  [Y (float)]  [Z (float)]</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Source number 100 with the coordinates (12.3 [m], 100.5 [m], 20.3 [m]).</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">100  12.3  100.5  20.3</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setText(QtGui.QApplication.translate("generate_survey", "Source\n"
"File [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_src.setText(QtGui.QApplication.translate("generate_survey", "Browse", None, QtGui.QApplication.UnicodeUTF8))

