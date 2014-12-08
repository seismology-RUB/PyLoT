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
        try:
            st = data.select(component=comp)
            numStations = len(st)
            for n in range(numStations):
                stat = st[n].stats.station
                stationButtons.append(QPushButton('%s'.format(stat)))
        except:
            for n in range(5):
                stationButtons.append(QPushButton('ST{0:02d}'.format(n+1)))
        for button in stationButtons:
            layout.addWidget(button)
        return layout