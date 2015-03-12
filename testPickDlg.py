#!/usr/bin/env python

import sys
from PySide.QtGui import QApplication
from obspy.core import read
from pylot.core.util.widgets import PickDlg

app = QApplication(sys.argv)

data = read()
win = PickDlg(data=data)
win.show()
app.exec_()
