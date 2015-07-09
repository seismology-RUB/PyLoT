import sys
from PySide.QtCore import QThread, Signal

class WorkerThread(QThread):
    message = Signal(str)

    def __init__(self, func, data, param):
        super(WorkerThread, self).__init__()
        self.func = func
        self.data = data
        self.param = param

    def run(self):
        sys.stdout = self

        self.func(self.data, self.param)

    def write(self, text):
        self.message.emit(text)
