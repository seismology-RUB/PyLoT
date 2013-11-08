# -*- coding: utf-8 -*-
#--------------------------------------------------------
# Purpose: Convience imports for PyLoT
#
'''
PyLoT - the Python picking and Localization Tool

The Python picking and Localisation Tool

This python library contains a graphical user interfaces for picking
seismic phases. This software needs ObsPy (http://github.com/obspy/obspy/wiki)
and the Qt4 libraries to be installed first.

PILOT has been developed in Mathworks' MatLab. In order to distribute
PILOT without facing portability problems, it has been decided to re-
develop the software package in Python. The great work of the ObsPy
group allows easy handling of a bunch of seismic data and PyLoT will
benefit a lot compared to the former MatLab version.

The development of PyLoT is part of the joint research project MAGS2.

:copyright:
    The PyLoT-Development Team
:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
'''

import os.path as osp
from obspy.core import read
