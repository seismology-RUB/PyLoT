import sys
from PySide.QtCore import QThread, Signal

class WorkerThread(QThread):
    message = Signal(str)

    def __init__(self, parent, func, data, param):
        super(WorkerThread, self).__init__()
        self.setParent(parent)
        self.func = func
        self.data = data
        self.param = param

    def run(self):
        sys.stdout = self

        picks = self.func(self.data, self.param)

        try:
            self.parent().addPicks(picks)
        except AttributeError:
            print picks

    def write(self, text):
        self.message.emit(text)
