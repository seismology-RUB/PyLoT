#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylot.core.io.phases import writephases
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def export(picks, fnout, parameter, eventinfo):
    '''
    Take <picks> dictionary and exports picking data to a HASH
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param fnout: complete path to the exporting obs file
    :type fnout: str
    
    :param: parameter, all input information
    :type:  object

    :param: eventinfo, source information needed for HASH format
    :type:  list object
    '''
    # write phases to HASH-phase file
    writephases(picks, 'HASH', fnout, parameter, eventinfo)
