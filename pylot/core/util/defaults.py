#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 12:31:25 2014

@author: sebastianw
"""

import os
import platform

from pylot.core.util.utils import readDefaultFilterInformation
from pylot.core.loc import hypo71
from pylot.core.loc import hypodd
from pylot.core.loc import hyposat
from pylot.core.loc import nll
from pylot.core.loc import velest


# determine system dependent path separator
system_name = platform.system()
if system_name in ["Linux", "Darwin"]:
    SEPARATOR = '/'
elif system_name == "Windows":
    SEPARATOR = '\\'

# suffix for phase name if not phase identified by last letter (P, p, etc.)
ALTSUFFIX = ['diff', 'n', 'g', '1', '2', '3']

FILTERDEFAULTS = readDefaultFilterInformation(os.path.join(os.path.expanduser('~'),
                                                           '.pylot',
                                                           'pylot.in'))

TIMEERROR_DEFAULTS = os.path.join(os.path.expanduser('~'),
                                  '.pylot',
                                  'PILOT_TimeErrors.in')

OUTPUTFORMATS = {'.xml': 'QUAKEML',
                 '.cnv': 'CNV',
                 '.obs': 'NLLOC_OBS'}

LOCTOOLS = dict(nll=nll, hyposat=hyposat, velest=velest, hypo71=hypo71, hypodd=hypodd)


class SetChannelComponents(object):
    def __init__(self):
        self.setDefaultCompPosition()

    def setDefaultCompPosition(self):
        # default component order
        self.compPosition_Map = dict(Z=2, N=1, E=0)
        self.compName_Map = {'3': 'Z',
                             '1': 'N',
                             '2': 'E'}

    def _getCurrentPosition(self, component):
        for key, value in self.compName_Map.items():
            if value == component:
                return key, value
        errMsg = 'getCurrentPosition: Could not find former position of component {}.'.format(component)
        raise ValueError(errMsg)

    def _switch(self, component, component_alter):
        # Without switching, multiple definitions of the same alter_comp are possible
        old_alter_comp, _ = self._getCurrentPosition(component)
        old_comp = self.compName_Map[component_alter]
        if not old_alter_comp == component_alter and not old_comp == component:
            self.compName_Map[old_alter_comp] = old_comp
            print('switch: Automatically switched component {} to {}'.format(old_alter_comp, old_comp))

    def setCompPosition(self, component_alter, component, switch=True):
        component_alter = str(component_alter)
        if not component_alter in self.compName_Map.keys():
            errMsg = 'setCompPosition: Unrecognized alternative component {}. Expecting one of {}.'
            raise ValueError(errMsg.format(component_alter, self.compName_Map.keys()))
        if not component in self.compPosition_Map.keys():
            errMsg = 'setCompPosition: Unrecognized target component {}. Expecting one of {}.'
            raise ValueError(errMsg.format(component, self.compPosition_Map.keys()))
        print('setCompPosition: set component {} to {}'.format(component_alter, component))
        if switch:
            self._switch(component, component_alter)
        self.compName_Map[component_alter] = component

    def getCompPosition(self, component):
        return self._getCurrentPosition(component)[0]

    def getPlotPosition(self, component):
        component = str(component)
        if component in self.compPosition_Map.keys():
            return self.compPosition_Map[component]
        elif component in self.compName_Map.keys():
            return self.compPosition_Map[self.compName_Map[component]]
        else:
            errMsg = 'getCompPosition: Unrecognized component {}. Expecting one of {} or {}.'
            raise ValueError(errMsg.format(component, self.compPosition_Map.keys(), self.compName_Map.keys()))
