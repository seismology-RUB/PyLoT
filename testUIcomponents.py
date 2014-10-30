#!/usr/bin/env python

import sys, time
from PySide.QtGui import QApplication
from pylot.core.util.widgets import FilterOptionsDialog, PropertiesDlg, HelpForm

dialogs = [FilterOptionsDialog, PropertiesDlg, HelpForm]

app = QApplication(sys.argv)

for dlg in dialogs:
	win = dlg()
	win.show()
	time.sleep(10)
	win.destroy()