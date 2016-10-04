#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylot.core.io.phases import writephases
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()

def export(picks, fnout):
    '''
    Take <picks> dictionary and exports picking data to a NLLOC-obs
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param fnout: complete path to the exporting obs file
    :type fnout: str
    '''
    # write phases to NLLoc-phase file
    writephases(picks, 'HYPO71', fnout)
