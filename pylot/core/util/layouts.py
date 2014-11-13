#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
'''
Created on 10.11.2014

@author: sebastianw
'''

from PySide.QtGui import (QVBoxLayout,
                          QPushButton)

def layoutStationButtons(data, comp):
        layout = QVBoxLayout()
        stationButtons = []
        st = data.select(component=comp)
        numStations = len(st)
        for n in range(numStations):
            stationButtons[n] = QPushButton('%s'.format(
                                            st[n].stats.station))
        layout.addWidget(stationButtons)

        return layout