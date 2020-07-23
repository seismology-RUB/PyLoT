#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 12:31:25 2014

@author: sebastianw
"""

import os
import platform

from pylot.core.loc import hypo71
from pylot.core.loc import hypodd
from pylot.core.loc import hyposat
from pylot.core.loc import nll
from pylot.core.loc import velest
from pylot.core.util.utils import readDefaultFilterInformation

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
                 '.obs': 'NLLOC_OBS',
                 '.focmec': 'FOCMEC'}

LOCTOOLS = dict(nll=nll, hyposat=hyposat, velest=velest, hypo71=hypo71, hypodd=hypodd)
