# Main program: QtPyLoT.py
#
# 
#
#
#
#
#
#
#
#
#
#

import os
import platform
import sys
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import helpform
from pylot.core.util import _getVersionString

# Version information
__version__ = _getVersionString()

class MainWindow(QMainWindow):
    
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        
        self.Stream = None


# Creating a Qt application
pylot_app = QApplication(sys.argv)

pylot_main = QWidget()
pylot_main.setWindowTitle('TestWindow')

ok_btn = QPushButton('OK', pylot_main)

@pyqtSlot()
def on_click():
	print('clicked')

@pyqtSlot()
def on_press():
	print('pressed')

@pyqtSlot()
def on_release():
	print('released')

# Connect signals to slots
ok_btn.clicked.connect(on_click)
ok_btn.pressed.connect(on_press)
ok_btn.pressed.connect(on_release)

# Show window and run the app
pylot_main.show()
pylot_app.exec_()
