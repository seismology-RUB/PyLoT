# -*- coding: utf-8 -*-
import sys
from PySide.QtCore import QThread, Signal

class AutoPickThread(QThread):
    message = Signal(str)

    def __init__(self, parent, func, data, param):
        super(AutoPickThread, self).__init__()
        self.setParent(parent)
        self.func = func
        self.data = data
        self.param = param

    def run(self):
        sys.stdout = self

        picks = self.func(self.data, self.param)

        print("Autopicking finished!\n")

        try:
            for station in picks:
                self.parent().addPicks(station, picks[station], type='auto')
        except AttributeError:
            print(picks)
        # plot picks to section
        self.parent().drawPicks(picktype='auto')


    def write(self, text):
        self.message.emit(text)
