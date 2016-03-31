#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy import read_events

from pylot.core.read.io import picks_from_evt
from pylot.core.util.pdf import ProbabilityDensityFunction
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'sebastianw'


def read_data(fn):
    """
    Reads pick data from QuakeML files named FN and returns a dictionary
    containing a ProbabilityDensityFunction object for each pick.
    :param fn: name of the QuakeML file which contains the picks
    :type fn: str
    :return: a dictionary containing the picks represented as pdfs
    """
    pdf_picks = picks_from_evt(read_events(fn)[0])

    for station, phases in pdf_picks.items():
        for phase, values in phases.items():
            phases[phase] = ProbabilityDensityFunction.fromPick(values['epp'],
                                                                values['mpp'],
                                                                values['lpp'],
                                                                type='exp')

    return pdf_picks
