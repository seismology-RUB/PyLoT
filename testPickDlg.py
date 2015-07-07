#!/usr/bin/env python

import sys
import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from PySide.QtGui import QApplication
from obspy.core import read
from pylot.core.util.widgets import PickDlg
import icons_rc

app = QApplication(sys.argv)

data = read()
win = PickDlg(data=data)
win.show()
app.exec_()
