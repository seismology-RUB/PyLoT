# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'generate_survey.ui'
#
# Created: Wed Jun 15 11:56:01 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_generate_survey(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(380, 160)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../asp3d_icon.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        Dialog.setWindowIcon(icon)
        self.verticalLayout = QtGui.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.gridLayout = QtGui.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.lineEdit_rec = QtGui.QLineEdit(Dialog)
        self.lineEdit_rec.setObjectName("lineEdit_rec")
        self.gridLayout.addWidget(self.lineEdit_rec, 0, 1, 1, 1)
        self.pushButton_rec = QtGui.QPushButton(Dialog)
        self.pushButton_rec.setObjectName("pushButton_rec")
        self.gridLayout.addWidget(self.pushButton_rec, 0, 2, 1, 1)
        self.label_rec = QtGui.QLabel(Dialog)
        self.label_rec.setObjectName("label_rec")
        self.gridLayout.addWidget(self.label_rec, 0, 0, 1, 1)
        self.lineEdit_obs = QtGui.QLineEdit(Dialog)
        self.lineEdit_obs.setObjectName("lineEdit_obs")
        self.gridLayout.addWidget(self.lineEdit_obs, 2, 1, 1, 1)
        self.label_obs = QtGui.QLabel(Dialog)
        self.label_obs.setObjectName("label_obs")
        self.gridLayout.addWidget(self.label_obs, 2, 0, 1, 1)
        self.pushButton_obs = QtGui.QPushButton(Dialog)
        self.pushButton_obs.setObjectName("pushButton_obs")
        self.gridLayout.addWidget(self.pushButton_obs, 2, 2, 1, 1)
        self.label_src = QtGui.QLabel(Dialog)
        self.label_src.setObjectName("label_src")
        self.gridLayout.addWidget(self.label_src, 1, 0, 1, 1)
        self.lineEdit_src = QtGui.QLineEdit(Dialog)
        self.lineEdit_src.setObjectName("lineEdit_src")
        self.gridLayout.addWidget(self.lineEdit_src, 1, 1, 1, 1)
        self.pushButton_src = QtGui.QPushButton(Dialog)
        self.pushButton_src.setObjectName("pushButton_src")
        self.gridLayout.addWidget(self.pushButton_src, 1, 2, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.buttonBox = QtGui.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), Dialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QtGui.QApplication.translate("Dialog", "Generate new Survey", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_rec.setText(QtGui.QApplication.translate("Dialog", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_rec.setText(QtGui.QApplication.translate("Dialog", "Receiver\n"
"File", None, QtGui.QApplication.UnicodeUTF8))
        self.label_obs.setText(QtGui.QApplication.translate("Dialog", "Seismogram\n"
"Directory", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_obs.setText(QtGui.QApplication.translate("Dialog", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setText(QtGui.QApplication.translate("Dialog", "Source\n"
"File", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_src.setText(QtGui.QApplication.translate("Dialog", "Browse", None, QtGui.QApplication.UnicodeUTF8))

