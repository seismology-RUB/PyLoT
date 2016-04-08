#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy

import numpy as np

import matplotlib.pyplot as plt

from obspy import read_events

from pylot.core.read.io import picks_from_evt
from pylot.core.util.pdf import ProbabilityDensityFunction
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'sebastianw'


class Comparison(object):
    """
    A Comparison object contains information on the evaluated picks' probability
    density function and compares these in terms of building the difference of
    compared pick sets. The results can be displayed as histograms showing its
    properties.
    """

    def __init__(self, **kwargs):
        names = list()
        self._pdfs = dict()
        for name, fn in kwargs.items():
            self._pdfs[name] = PDFDictionary.from_quakeml(fn)
            names.append(name)
        if len(names) > 2:
            raise ValueError('Comparison is only defined for two '
                             'arguments!')
        self._names = names
        self._compare = self.compare_picksets()

    def __nonzero__(self):
        if not len(self.names) == 2 or not self._pdfs:
            return False
        return True

    def get(self, name):
        return self._pdfs[name]

    @property
    def names(self):
        return self._names

    @names.setter
    def names(self, names):
        assert isinstance(names, list) and len(names) == 2, 'variable "names"' \
                                                            ' is either not a' \
                                                            ' list or its ' \
                                                            'length is not 2:' \
                                                            'names : {names}'.format(
            names=names)
        self._names = names

    @property
    def comparison(self):
        return self._compare

    @property
    def stations(self):
        return self.comparison.keys()

    @property
    def nstations(self):
        return len(self.stations)

    def compare_picksets(self, type='exp'):
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
        compare_pdfs = dict()

        pdf_a = self.get(self.names[0]).pdf_data(type)
        pdf_b = self.get(self.names[1]).pdf_data(type)

        for station, phases in pdf_a.items():
            if station in pdf_b.keys():
                compare_pdf = dict()
                for phase in phases:
                    if phase in pdf_b[station].keys():
                        compare_pdf[phase] = phases[phase] - pdf_b[station][
                            phase]
                if compare_pdf is not None:
                    compare_pdfs[station] = compare_pdf

        return compare_pdfs

    def plot(self):
        nstations = self.nstations
        stations = self.stations
        istations = range(nstations)
        fig, axarr = plt.subplots(nstations, 2, sharex='col', sharey='row')

        for n in istations:
            station = stations[n]
            compare_pdf = self.comparison[station]
            for l, phase in enumerate(compare_pdf.keys()):
                axarr[n, l].plot(compare_pdf[phase].axis,
                                 compare_pdf[phase].data)
                if n is 0:
                    axarr[n, l].set_title(phase)
                if l is 0:
                    axann = axarr[n, l].annotate(station, xy=(.05, .5),
                                                 xycoords='axes fraction')
                    bbox_props = dict(boxstyle='round', facecolor='lightgrey',
                                      alpha=.7)
                    axann.set_bbox(bbox_props)
                if n == int(np.median(istations)) and l is 0:
                    label = 'probability density (qualitative)'
                    axarr[n, l].set_ylabel(label)
        plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
        plt.setp([a.get_yticklabels() for a in axarr[:, 0]], visible=False)

        plt.show()

    def get_expectation_array(self, phase):
        pdf_dict = self.comparison
        exp_array = list()
        for station, phases in pdf_dict.items():
            exp_array.append(phases[phase].expectation())
        return exp_array

    def get_std_array(self, phase):
        pdf_dict = self.comparison
        std_array = list()
        for station, phases in pdf_dict.items():
            std_array.append(phases[phase].standard_deviation())
        return std_array

    def hist_expectation(self, phases='all', bins=20, normed=False):
        phases.strip()
        if phases.find('all') is 0:
            phases = 'ps'
        phases = phases.upper()
        nsp = len(phases)
        fig, axarray = plt.subplots(1, nsp, sharey=True)
        for n, phase in enumerate(phases):
            ax = axarray[n]
            data = self.get_expectation_array(phase)
            xlims = [min(data), max(data)]
            ax.hist(data, range=xlims, bins=bins, normed=normed)
            title_str = 'phase: {0}, samples: {1}'.format(phase, len(data))
            ax.set_title(title_str)
            ax.set_xlabel('expectation [s]')
            if n is 0:
                ax.set_ylabel('abundance [-]')
        plt.setp([a.get_yticklabels() for a in axarray[1:]], visible=False)
        plt.show()

    def hist_standard_deviation(self, phases='all', bins=20, normed=False):
        phases.strip()
        if phases.find('all') == 0:
            phases = 'ps'
        phases = phases.upper()
        nsp = len(phases)
        fig, axarray = plt.subplots(1, nsp, sharey=True)
        for n, phase in enumerate(phases):
            ax = axarray[n]
            data = self.get_std_array(phase)
            xlims = [min(data), max(data)]
            ax.hist(data, range=xlims, bins=bins, normed=normed)
            title_str = 'phase: {0}, samples: {1}'.format(phase, len(data))
            ax.set_title(title_str)
            ax.set_xlabel('standard deviation [s]')
            if n is 0:
                ax.set_ylabel('abundance [-]')
        plt.setp([a.get_yticklabels() for a in axarray[1:]], visible=False)
        plt.show()

    def hist(self, type='std'):
        pass


class PDFDictionary(object):
    """
    A PDFDictionary is a dictionary like object containing structured data on
    the probability density function of seismic phase onsets.
    """

    def __init__(self, data):
        self._pickdata = data

    def __nonzero__(self):
        if len(self.pick_data) < 1:
            return False
        else:
            return True

    @property
    def pick_data(self):
        return self._pickdata

    @pick_data.setter
    def pick_data(self, data):
        self._pickdata = data

    @classmethod
    def from_quakeml(self, fn):
        cat = read_events(fn)
        if len(cat) > 1:
            raise NotImplementedError('reading more than one event at the same '
                                      'time is not implemented yet! Sorry!')
        return PDFDictionary(picks_from_evt(cat[0]))

    def pdf_data(self, type='exp'):
        """
        Returns probabiliy density function dictionary containing the
        representation of the actual pick_data.
        :param type: type of the returned
         `~pylot.core.util.pdf.ProbabilityDensityFunction` object
        :type type: str
        :return: a dictionary containing the picks represented as pdfs
        """

        pdf_picks = copy.deepcopy(self.pick_data)

        for station, phases in pdf_picks.items():
            for phase, values in phases.items():
                phases[phase] = ProbabilityDensityFunction.fromPick(
                    values['epp'],
                    values['mpp'],
                    values['lpp'],
                    type=type)

        return pdf_picks


comp_obj = Comparison(manual='/home/sebastianp/Data/Insheim/e0057.026.13.xml', auto='/home/sebastianp/Data/Insheim/e0057.026.13.xml')
comp_obj.plot()
comp_obj.hist_expectation()
comp_obj.hist_standard_deviation()