#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 17:47:25 2015

@author: sebastianw
"""

from pylot.core.io.data import SeiscompDataStructure, PilotDataStructure

DATASTRUCTURE = {'PILOT': PilotDataStructure, 'SeisComP': SeiscompDataStructure,
                 None: None}
