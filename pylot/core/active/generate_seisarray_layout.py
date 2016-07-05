# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'generate_seisarray_layout.ui'
#
# Created: Tue Jul  5 13:55:55 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_generate_seisarray(object):
    def setupUi(self, generate_seisarray):
        generate_seisarray.setObjectName("generate_seisarray")
        generate_seisarray.resize(403, 221)
        generate_seisarray.setMinimumSize(QtCore.QSize(400, 220))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("../asp3d_icon.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        generate_seisarray.setWindowIcon(icon)
        self.verticalLayout_5 = QtGui.QVBoxLayout(generate_seisarray)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_rec = QtGui.QLabel(generate_seisarray)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_rec.sizePolicy().hasHeightForWidth())
        self.label_rec.setSizePolicy(sizePolicy)
        self.label_rec.setToolTip("")
        self.label_rec.setObjectName("label_rec")
        self.horizontalLayout.addWidget(self.label_rec)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.radioButton_normal = QtGui.QRadioButton(generate_seisarray)
        self.radioButton_normal.setChecked(True)
        self.radioButton_normal.setObjectName("radioButton_normal")
        self.horizontalLayout.addWidget(self.radioButton_normal)
        self.radioButton_interpolatable = QtGui.QRadioButton(generate_seisarray)
        self.radioButton_interpolatable.setObjectName("radioButton_interpolatable")
        self.horizontalLayout.addWidget(self.radioButton_interpolatable)
        self.verticalLayout_4.addLayout(self.horizontalLayout)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.lineEdit_rec = QtGui.QLineEdit(generate_seisarray)
        self.lineEdit_rec.setObjectName("lineEdit_rec")
        self.horizontalLayout_4.addWidget(self.lineEdit_rec)
        self.pushButton_rec = QtGui.QPushButton(generate_seisarray)
        self.pushButton_rec.setObjectName("pushButton_rec")
        self.horizontalLayout_4.addWidget(self.pushButton_rec)
        self.verticalLayout_4.addLayout(self.horizontalLayout_4)
        self.verticalLayout.addLayout(self.verticalLayout_4)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_src = QtGui.QLabel(generate_seisarray)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_src.sizePolicy().hasHeightForWidth())
        self.label_src.setSizePolicy(sizePolicy)
        self.label_src.setObjectName("label_src")
        self.verticalLayout_2.addWidget(self.label_src)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.lineEdit_src = QtGui.QLineEdit(generate_seisarray)
        self.lineEdit_src.setObjectName("lineEdit_src")
        self.horizontalLayout_5.addWidget(self.lineEdit_src)
        self.pushButton_src = QtGui.QPushButton(generate_seisarray)
        self.pushButton_src.setObjectName("pushButton_src")
        self.horizontalLayout_5.addWidget(self.pushButton_src)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.verticalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_obs = QtGui.QLabel(generate_seisarray)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_obs.sizePolicy().hasHeightForWidth())
        self.label_obs.setSizePolicy(sizePolicy)
        self.label_obs.setObjectName("label_obs")
        self.verticalLayout_3.addWidget(self.label_obs)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.lineEdit_pts = QtGui.QLineEdit(generate_seisarray)
        self.lineEdit_pts.setObjectName("lineEdit_pts")
        self.horizontalLayout_6.addWidget(self.lineEdit_pts)
        self.pushButton_obs = QtGui.QPushButton(generate_seisarray)
        self.pushButton_obs.setObjectName("pushButton_obs")
        self.horizontalLayout_6.addWidget(self.pushButton_obs)
        self.verticalLayout_3.addLayout(self.horizontalLayout_6)
        self.verticalLayout.addLayout(self.verticalLayout_3)
        self.verticalLayout_5.addLayout(self.verticalLayout)
        self.buttonBox = QtGui.QDialogButtonBox(generate_seisarray)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_5.addWidget(self.buttonBox)

        self.retranslateUi(generate_seisarray)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), generate_seisarray.reject)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), generate_seisarray.accept)
        QtCore.QMetaObject.connectSlotsByName(generate_seisarray)

    def retranslateUi(self, generate_seisarray):
        generate_seisarray.setWindowTitle(QtGui.QApplication.translate("generate_seisarray", "Generate new Seismic Array", None, QtGui.QApplication.UnicodeUTF8))
        self.label_rec.setText(QtGui.QApplication.translate("generate_seisarray", "Measured Receivers", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_normal.setToolTip(QtGui.QApplication.translate("generate_seisarray", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Load <span style=\" font-weight:600;\">normal </span>measured receiver input file. The input file must be in the following format:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Containing in each line, separated by spaces:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">[trace ID (int)]  [X (float)]  [Y (float)]  [Z (float)]</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Geophone with the trace ID 50 and the coordinates (10.5 [m], 20.4 [m], 30.3 [m]).</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">50  10.5  20.4  30.3</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_normal.setText(QtGui.QApplication.translate("generate_seisarray", "normal [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_interpolatable.setToolTip(QtGui.QApplication.translate("generate_seisarray", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Load<span style=\" font-weight:600;\"> </span>measured receiver input file that can be <span style=\" font-weight:600;\">interpolated</span>. The input file must be in the following format:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Containing in each line, separated by spaces:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">[trace ID (int)] [receiver line ID (int)] [number of the geophone on receiver line (int)]  [X (float)]  [Y (float)]  [Z (float)]</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Third geophone on the second receiver line with the trace ID 50 and the coordinates (10.5 [m], 20.4 [m], 30.3 [m]).</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">50  2  3  10.5  20.4  30.3</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_interpolatable.setText(QtGui.QApplication.translate("generate_seisarray", "interpolatable [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_rec.setText(QtGui.QApplication.translate("generate_seisarray", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setToolTip(QtGui.QApplication.translate("generate_seisarray", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Load measured sources input file to improve interpolation precision. The input file must be in the following format:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Containing in each line, separated by spaces:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">[source ID (int)] [X (float)]  [Y (float)]  [Z (float)]</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Shot number 100 with the coordinates (12.3 [m], 100.5 [m], 20.3 [m]).</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">100  12.3  100.5  20.3</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setText(QtGui.QApplication.translate("generate_seisarray", "Measured Sources (optional) [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_src.setText(QtGui.QApplication.translate("generate_seisarray", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_obs.setToolTip(QtGui.QApplication.translate("generate_seisarray", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Load measured points input file to improve interpolation precision. The input file must be in the following format:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Containing in each line, separated by spaces:</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">[point ID (int)] [X (float)]  [Y (float)]  [Z (float)]</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-weight:600;\">For example:</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Point number 100 with the coordinates (12.3 [m], 100.5 [m], 20.3 [m]).</p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">100  12.3  100.5  20.3</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.label_obs.setText(QtGui.QApplication.translate("generate_seisarray", "Additional measured points (optional) [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_obs.setText(QtGui.QApplication.translate("generate_seisarray", "Browse", None, QtGui.QApplication.UnicodeUTF8))

