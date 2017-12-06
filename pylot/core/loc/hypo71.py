#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylot.core.io.phases import writephases
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def export(picks, fnout, parameter):
    '''
    Take <picks> dictionary and exports picking data to a HYPO71
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param fnout: complete path to the exporting obs file
    :type fnout: str

    :param parameter: all input information
    :type parameter:  object
    '''
    # write phases to HYPO71-phase file
    writephases(picks, 'HYPO71', fnout, parameter)
