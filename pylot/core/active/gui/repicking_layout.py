# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'repicking_layout.ui'
#
# Created: Mon Aug 29 10:26:23 2016
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_repicking(object):
    def setupUi(self, repicking):
        repicking.setObjectName("repicking")
        repicking.resize(640, 480)
        self.verticalLayout = QtGui.QVBoxLayout(repicking)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pushButton_repick = QtGui.QPushButton(repicking)
        self.pushButton_repick.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_repick.setObjectName("pushButton_repick")
        self.horizontalLayout.addWidget(self.pushButton_repick)
        self.pushButton_delete = QtGui.QPushButton(repicking)
        self.pushButton_delete.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_delete.setObjectName("pushButton_delete")
        self.horizontalLayout.addWidget(self.pushButton_delete)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout_plot = QtGui.QVBoxLayout()
        self.verticalLayout_plot.setObjectName("verticalLayout_plot")
        self.verticalLayout.addLayout(self.verticalLayout_plot)
        self.buttonBox = QtGui.QDialogButtonBox(repicking)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Close)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(repicking)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), repicking.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), repicking.reject)
        QtCore.QMetaObject.connectSlotsByName(repicking)

    def retranslateUi(self, repicking):
        repicking.setWindowTitle(QtGui.QApplication.translate("repicking", "Repicking", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_repick.setText(QtGui.QApplication.translate("repicking", "Repick", None, QtGui.QApplication.UnicodeUTF8))
        self.pushButton_delete.setText(QtGui.QApplication.translate("repicking", "Delete", None, QtGui.QApplication.UnicodeUTF8))

