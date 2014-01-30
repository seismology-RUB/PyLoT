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
        
        filterDockWidget = QDockWidget("Filter Options", self)
    
 
 class PickWindow(QDialog):

	def __init__(self, station=None, parent=None):
		super(PickWindow, self).__init__(parent)
		
		filterDockWidget = FilterOptionsDock()

class PropertiesWindow(QDialog):
		
	def __init__(self, parent=None):
		super(PropertiesWindow, self).__init__(parent)

class FilterOptionsDock(QDockWidget):
	
	def __init__(self, filterOptions=None):
		super(FilterOptionsDock, self).__init__()
		
		if filterOptions and not isinstance(filterOptions, FilterOptions):
			try:
                     fOptions = FilterOptions(filterOptions)
                 except e:
                     raise OptionsError, '%s' % e

class FilterOptions(object):
	
	def __init__(self, filtertype=None, freq=None, order=None):
         self.__filterInformation = {}
         self._setfilterType(filtertype)
         self._setFreq(freq)
         self._setOrder(order)
     
     def _getFreq(self):
         return self.__filterInformation['freq']
         
     def _setFreq(self, freq):
         self.__filterInformation['freq'] = freq
     
     def _getOrder(self):
         return self.__filterInformation['order']
         
     def _setOrder(self, order):
         self.__filterInformation['order'] = order
         
     def _getFilterType(self):
         return self.__filterInformation['filtertype']
     
     def _setFilterType(self, filtertype):
         self.__filterInformation['filtertype'] = filtertype         
         
     filterType = property(fget=_getFilterType, fset=_setFilterType)
     order = property(fget=_getOrder, fset=_setOrder)
     freq = property(fget=_getFreq, fset=_setFreq)
     
class OptionsError(Exception): pass        

if __name__ == '__main__':
    ##Creating a Qt application
    pylot_app = QApplication(sys.argv)

    pylot_main = MainWindow()
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
