#!/usr/bin/env python

import sys, time
from PySide.QtGui import QApplication
from pylot.core.util.widgets import PropertiesDlg

app = QApplication(sys.argv)

win = PropertiesDlg()
win.show()
app.exec_()