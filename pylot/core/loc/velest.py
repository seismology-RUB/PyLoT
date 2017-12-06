#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylot.core.io.phases import writephases
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def export(picks, fnout, eventinfo, parameter=None):
    '''
    Take <picks> dictionary and exports picking data to a VELEST-cnv
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param fnout: complete path to the exporting obs file
    :type fnout: str
    
    :param eventinfo: source time needed for VELEST-cnv format
    :type eventinfo:  list object

    :param parameter: all input information
    :type parameter:  object
    '''
    # write phases to VELEST-phase file
    writephases(picks, 'VELEST', fnout, parameter, eventinfo)
