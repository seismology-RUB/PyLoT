#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy import read_events

from pylot.core.read.io import picks_from_evt
from pylot.core.util.pdf import ProbabilityDensityFunction
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'sebastianw'


def read_data(fn, type='exp'):
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
                                                                type=type)

    return pdf_picks


def compare_picksets(a, b):
    """
    Compare two picksets A and B and return a dictionary compiling the results.
    Comparison is carried out with the help of pdf representation of the picks
    and a probabilistic approach to the time difference of two onset
    measurements.
    :param a: filename for pickset A
    :type a: str
    :param b: filename for pickset B
    :type b: str
    :return: dictionary containing the resulting comparison pdfs for all picks
    :rtype: dict
    """
    pdf_a = read_data(a)
    pdf_b = read_data(b)

    compare_pdfs = dict()

    for station, phases in pdf_a.items():
        if station in pdf_b.keys():
            compare_pdf = dict()
            for phase in phases:
                if phase in pdf_b[station].keys():
                    compare_pdf[phase] = phases[phase] - pdf_b[station][phase]
            if compare_pdf is not None:
                compare_pdfs[station] = compare_pdf

    return compare_pdfs
