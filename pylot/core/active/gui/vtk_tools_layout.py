# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'vtk_tools_layout.ui'
#
# Created: Thu Aug  4 12:07:57 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_vtk_tools(object):
    def setupUi(self, vtk_tools):
        vtk_tools.setObjectName("vtk_tools")
        vtk_tools.resize(422, 471)
        self.verticalLayout_9 = QtGui.QVBoxLayout(vtk_tools)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.verticalLayout = QtGui.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.verticalLayout_4 = QtGui.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_rec = QtGui.QLabel(vtk_tools)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_rec.sizePolicy().hasHeightForWidth())
        self.label_rec.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_rec.setFont(font)
        self.label_rec.setToolTip("")
        self.label_rec.setObjectName("label_rec")
        self.horizontalLayout.addWidget(self.label_rec)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.radioButton_abs = QtGui.QRadioButton(vtk_tools)
        self.radioButton_abs.setChecked(True)
        self.radioButton_abs.setObjectName("radioButton_abs")
        self.horizontalLayout.addWidget(self.radioButton_abs)
        self.radioButton_rel = QtGui.QRadioButton(vtk_tools)
        self.radioButton_rel.setObjectName("radioButton_rel")
        self.horizontalLayout.addWidget(self.radioButton_rel)
        self.verticalLayout_4.addLayout(self.horizontalLayout)
        self.verticalLayout_6 = QtGui.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.label_3 = QtGui.QLabel(vtk_tools)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_6.addWidget(self.label_3)
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.lineEdit_vg = QtGui.QLineEdit(vtk_tools)
        self.lineEdit_vg.setObjectName("lineEdit_vg")
        self.horizontalLayout_2.addWidget(self.lineEdit_vg)
        self.pushButton_vg = QtGui.QPushButton(vtk_tools)
        self.pushButton_vg.setObjectName("pushButton_vg")
        self.horizontalLayout_2.addWidget(self.pushButton_vg)
        self.verticalLayout_6.addLayout(self.horizontalLayout_2)
        self.verticalLayout_4.addLayout(self.verticalLayout_6)
        self.verticalLayout_7 = QtGui.QVBoxLayout()
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.label_4 = QtGui.QLabel(vtk_tools)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_7.addWidget(self.label_4)
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.lineEdit_vgref = QtGui.QLineEdit(vtk_tools)
        self.lineEdit_vgref.setEnabled(False)
        self.lineEdit_vgref.setObjectName("lineEdit_vgref")
        self.horizontalLayout_4.addWidget(self.lineEdit_vgref)
        self.pushButton_vgref = QtGui.QPushButton(vtk_tools)
        self.pushButton_vgref.setEnabled(False)
        self.pushButton_vgref.setObjectName("pushButton_vgref")
        self.horizontalLayout_4.addWidget(self.pushButton_vgref)
        self.verticalLayout_7.addLayout(self.horizontalLayout_4)
        self.verticalLayout_8 = QtGui.QVBoxLayout()
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.verticalLayout_7.addLayout(self.verticalLayout_8)
        self.label_5 = QtGui.QLabel(vtk_tools)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_7.addWidget(self.label_5)
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.lineEdit_vgout = QtGui.QLineEdit(vtk_tools)
        self.lineEdit_vgout.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.lineEdit_vgout.setObjectName("lineEdit_vgout")
        self.horizontalLayout_3.addWidget(self.lineEdit_vgout)
        self.pushButton_vtkout = QtGui.QPushButton(vtk_tools)
        self.pushButton_vtkout.setObjectName("pushButton_vtkout")
        self.horizontalLayout_3.addWidget(self.pushButton_vtkout)
        self.pushButton_parav = QtGui.QPushButton(vtk_tools)
        self.pushButton_parav.setObjectName("pushButton_parav")
        self.horizontalLayout_3.addWidget(self.pushButton_parav)
        self.verticalLayout_7.addLayout(self.horizontalLayout_3)
        self.start_vg = QtGui.QPushButton(vtk_tools)
        self.start_vg.setEnabled(False)
        self.start_vg.setObjectName("start_vg")
        self.verticalLayout_7.addWidget(self.start_vg)
        self.verticalLayout_4.addLayout(self.verticalLayout_7)
        self.verticalLayout.addLayout(self.verticalLayout_4)
        self.verticalLayout_2 = QtGui.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.line = QtGui.QFrame(vtk_tools)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout_2.addWidget(self.line)
        self.label_src = QtGui.QLabel(vtk_tools)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_src.sizePolicy().hasHeightForWidth())
        self.label_src.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_src.setFont(font)
        self.label_src.setObjectName("label_src")
        self.verticalLayout_2.addWidget(self.label_src)
        self.verticalLayout_3 = QtGui.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label = QtGui.QLabel(vtk_tools)
        self.label.setObjectName("label")
        self.verticalLayout_3.addWidget(self.label)
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.lineEdit_rays = QtGui.QLineEdit(vtk_tools)
        self.lineEdit_rays.setObjectName("lineEdit_rays")
        self.horizontalLayout_6.addWidget(self.lineEdit_rays)
        self.pushButton_rays = QtGui.QPushButton(vtk_tools)
        self.pushButton_rays.setObjectName("pushButton_rays")
        self.horizontalLayout_6.addWidget(self.pushButton_rays)
        self.verticalLayout_3.addLayout(self.horizontalLayout_6)
        self.verticalLayout_2.addLayout(self.verticalLayout_3)
        self.verticalLayout_5 = QtGui.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.label_2 = QtGui.QLabel(vtk_tools)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_5.addWidget(self.label_2)
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.lineEdit_raysout = QtGui.QLineEdit(vtk_tools)
        self.lineEdit_raysout.setObjectName("lineEdit_raysout")
        self.horizontalLayout_5.addWidget(self.lineEdit_raysout)
        self.pushButton_raysout = QtGui.QPushButton(vtk_tools)
        self.pushButton_raysout.setObjectName("pushButton_raysout")
        self.horizontalLayout_5.addWidget(self.pushButton_raysout)
        self.verticalLayout_5.addLayout(self.horizontalLayout_5)
        self.start_rays = QtGui.QPushButton(vtk_tools)
        self.start_rays.setEnabled(False)
        self.start_rays.setObjectName("start_rays")
        self.verticalLayout_5.addWidget(self.start_rays)
        self.verticalLayout_2.addLayout(self.verticalLayout_5)
        self.verticalLayout.addLayout(self.verticalLayout_2)
        self.verticalLayout_9.addLayout(self.verticalLayout)
        self.buttonBox = QtGui.QDialogButtonBox(vtk_tools)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_9.addWidget(self.buttonBox)

        self.retranslateUi(vtk_tools)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), vtk_tools.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), vtk_tools.reject)
        QtCore.QMetaObject.connectSlotsByName(vtk_tools)

    def retranslateUi(self, vtk_tools):
        vtk_tools.setWindowTitle(QtGui.QApplication.translate("vtk_tools", "VTK tools", None, QtGui.QApplication.UnicodeUTF8))
        self.label_rec.setText(QtGui.QApplication.translate("vtk_tools", "Velocity grid to VTK", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_abs.setToolTip(QtGui.QApplication.translate("vtk_tools", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
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
        self.radioButton_abs.setText(QtGui.QApplication.translate("vtk_tools", "absolute [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_rel.setToolTip(QtGui.QApplication.translate("vtk_tools", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
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
        self.radioButton_rel.setText(QtGui.QApplication.translate("vtk_tools", "relative [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("vtk_tools", "Browse for velocity grid file (\'vgrids.in\'):", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_vg.setText(QtGui.QApplication.translate("vtk_tools", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("vtk_tools", "Browse for reference velocity grid file (\'vgridsref.in\'):", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_vgref.setText(QtGui.QApplication.translate("vtk_tools", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("vtk_tools", "Output Filename:", None, QtGui.QApplication.UnicodeUTF8))
        self.lineEdit_vgout.setText(QtGui.QApplication.translate("vtk_tools", "vgrids.vtk", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_vtkout.setText(QtGui.QApplication.translate("vtk_tools", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_parav.setToolTip(QtGui.QApplication.translate("vtk_tools", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'Sans\'; font-size:10pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\">Call Paraview (Shell command: <span style=\" font-style:italic;\">paraview</span>) for the specified output filename in the current directory.</p></body></html>", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_parav.setText(QtGui.QApplication.translate("vtk_tools", "<- Paraview", None, QtGui.QApplication.UnicodeUTF8))
        self.start_vg.setText(QtGui.QApplication.translate("vtk_tools", "Start", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setToolTip(QtGui.QApplication.translate("vtk_tools", "Create VTK files from the FMTOMO output file \'rays.dat\'.\n"
"This will generate one file of ray paths for each individual source.", None, QtGui.QApplication.UnicodeUTF8))
        self.label_src.setText(QtGui.QApplication.translate("vtk_tools", "Rays to VTK [?]", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("vtk_tools", "Browse for input file (\'rays.dat\'):", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_rays.setText(QtGui.QApplication.translate("vtk_tools", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("vtk_tools", "Specify output directory:", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_raysout.setText(QtGui.QApplication.translate("vtk_tools", "Browse", None, QtGui.QApplication.UnicodeUTF8))
        self.start_rays.setText(QtGui.QApplication.translate("vtk_tools", "Start", None, QtGui.QApplication.UnicodeUTF8))

