# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 18:40:55 2014

@author: sebastianw
"""

from obspy.core.event import Pick


class ReferencePick(Pick):

    def __init__(self):
        Pick.__init__()
