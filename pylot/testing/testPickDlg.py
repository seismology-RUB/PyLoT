#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import matplotlib

matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4'] = 'PySide'

from PySide.QtGui import QApplication
from obspy.core import read
from pylot.core.util.widgets import PickDlg

app = QApplication(sys.argv)

data = read()
win = PickDlg(data=data)
win.show()
app.exec_()
