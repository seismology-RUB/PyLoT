# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 15:12:43 2013

@author: sebastianw
"""

from obspy.core import Trace as Obspytrace
from obspy.core.util import AttribDict


class Stats(AttribDict):
    pass


class Trace(Obspytrace):
    pass
         