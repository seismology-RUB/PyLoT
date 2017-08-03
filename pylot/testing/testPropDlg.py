#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from PySide.QtGui import QApplication
from pylot.core.util.widgets import PropertiesDlg

app = QApplication(sys.argv)

win = PropertiesDlg()
win.show()
app.exec_()
