#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylot.core.io.phases import writephases
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def export(picks, fnout, parameter):
    '''
    Take <picks> dictionary and exports picking data to a HYPOSAT
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param fnout: complete path to the exporting obs file
    :type fnout: str
    
    :param: parameter, all input information
    :type:  object
    '''
    # write phases to HYPOSAT-phase file
    writephases(picks, 'HYPOSAT', fnout, parameter)
