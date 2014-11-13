#!/usr/bin/env python

import sys, time
from PySide.QtGui import QApplication
from pylot.core.util.widgets import HelpForm

app = QApplication(sys.argv)

win = HelpForm()
win.show()
app.exec_()
