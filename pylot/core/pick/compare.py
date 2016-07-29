#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import os
import numpy as np
import glob
import matplotlib.pyplot as plt

from obspy import read_events

from pylot.core.io.phases import picksdict_from_picks
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
            if isinstance(fn, PDFDictionary):
                self._pdfs[name] = fn
            elif isinstance(fn, dict):
                self._pdfs[name] = PDFDictionary(fn)
            else:
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

        pdf_a = self.get(self.names[0]).generate_pdf_data(type)
        pdf_b = self.get(self.names[1]).generate_pdf_data(type)

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

    def plot(self, stations=None):
        if stations is None:
            nstations = self.nstations
            stations = self.stations
        else:
            nstations = len(stations)
        istations = range(nstations)
        fig, axarr = plt.subplots(nstations, 2, sharex='col', sharey='row')

        for n in istations:
            station = stations[n]
            if station not in self.comparison.keys():
                continue
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

    def get_all(self, phasename):
        pdf_dict = self.comparison
        rlist = list()
        for phases in pdf_dict.values():
            try:
                rlist.append(phases[phasename])
            except KeyError:
                continue
        return rlist

    def get_array(self, phase, method_name):
        import operator
        method = operator.methodcaller(method_name)
        pdf_list = self.get_all(phase)
        rarray = map(method, pdf_list)
        return np.array(rarray)

    def get_expectation_array(self, phase):
        return self.get_array(phase, 'expectation')

    def get_std_array(self, phase):
        return self.get_array(phase, 'standard_deviation')

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
        self._pdfdata = self.generate_pdf_data()

    def __nonzero__(self):
        if len(self.pick_data) < 1:
            return False
        else:
            return True

    def __getitem__(self, item):
        return self.pdf_data[item]

    @property
    def pdf_data(self):
        return self._pdfdata

    @pdf_data.setter
    def pdf_data(self, data):
        self._pdfdata = data

    @property
    def pick_data(self):
        return self._pickdata

    @pick_data.setter
    def pick_data(self, data):
        self._pickdata = data

    @property
    def stations(self):
        return self.pick_data.keys()

    @property
    def nstations(self):
        return len(self.stations)

    @classmethod
    def from_quakeml(self, fn):
        cat = read_events(fn)
        if len(cat) > 1:
            raise NotImplementedError('reading more than one event at the same '
                                      'time is not implemented yet! Sorry!')
        return PDFDictionary(picksdict_from_picks(cat[0]))

    def generate_pdf_data(self, type='exp'):
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
                if phase not in 'PS':
                    continue
                phases[phase] = ProbabilityDensityFunction.from_pick(
                    values['epp'],
                    values['mpp'],
                    values['lpp'],
                    type=type)

        return pdf_picks

    def plot(self, stations=None):
        '''
        plots the all probability density function for either desired STATIONS
        or all available date
        :param stations: list of stations to be plotted
        :type stations: list
        :return: matplotlib figure object containing the plot
        '''
        assert stations is not None or not isinstance(stations, list), \
            'parameter stations should be a list not {0}'.format(type(stations))
        if not stations:
            nstations = self.nstations
            stations = self.stations
        else:
            nstations = len(stations)

        istations = range(nstations)
        fig, axarr = plt.subplots(nstations, 2, sharex='col', sharey='row')
        hide_labels = True

        for n in istations:
            station = stations[n]
            pdfs = self.pdf_data[station]
            for l, phase in enumerate(pdfs.keys()):
                try:
                    axarr[n, l].plot(pdfs[phase].axis, pdfs[phase].data())
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
                except IndexError as e:
                    print('trying aligned plotting\n{0}'.format(e))
                    hide_labels = False
                    axarr[l].plot(pdfs[phase].axis, pdfs[phase].data())
                    axarr[l].set_title(phase)
                    if l is 0:
                        axann = axarr[l].annotate(station, xy=(.05, .5),
                                                     xycoords='axes fraction')
                        bbox_props = dict(boxstyle='round', facecolor='lightgrey',
                                          alpha=.7)
                        axann.set_bbox(bbox_props)
        if hide_labels:
            plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
            plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
            plt.setp([a.get_yticklabels() for a in axarr[:, 0]], visible=False)

        return fig


class PDFstatistics(object):
    '''
    To do:
    plots for std, quantiles,
    '''
    def __init__(self, directory):
        self.directory = directory
        self.arraylen = 0
        self.stations = {}
        self.p_std = {}
        self.s_std = {}
        self.theta0 = []
        self.theta1 = []
        self.theta2 = []


    def readTheta(self, arname, dir, fnpattern):
        exec('self.' + arname +' = []')
        filelist = glob.glob1(dir, fnpattern)
        for file in filelist:
            fid = open(os.path.join(dir,file), 'r')
            list = []
            for line in fid.readlines():
                list.append(eval(line))
            exec('self.' + arname + ' += list')
            fid.close()


    def fromFileList(self):
        self.makeFileList()
        self.getData()
        self.getStatistics()


    def makeFileList(self, fn_pattern='*'):
        self.evtlist = glob.glob1((os.path.join(self.directory)), '*.xml')


    def getData(self):
        for evt in self.evtlist:
            print evt
            self.stations[evt] = []
            self.p_std[evt] = []
            self.s_std[evt] = []
            self.getPDFDict(self.directory, evt)
            for station, pdfs in self.pdfdict.pdf_data.items():
                # print station, pdfs
                try:
                    p_std = pdfs['P'].standard_deviation()
                    self.theta0.append(pdfs['P'].qtile_dist_quot(0.015))
                    self.theta1.append(pdfs['P'].qtile_dist_quot(0.1))
                    self.theta2.append(pdfs['P'].qtile_dist_quot(0.2))
                except KeyError:
                    p_std = np.nan
                try:
                    s_std = pdfs['S'].standard_deviation()
                    self.theta0.append(pdfs['S'].qtile_dist_quot(0.015))
                    self.theta1.append(pdfs['S'].qtile_dist_quot(0.1))
                    self.theta2.append(pdfs['S'].qtile_dist_quot(0.2))
                except KeyError:
                    s_std = np.nan
                self.stations[evt].append(station)
                self.p_std[evt].append(p_std)
                self.s_std[evt].append(s_std)
                self.arraylen += 1
        self.makeArray()


    def makeArray(self):
        index = 0
        self.p_stdarray = np.zeros(self.arraylen)
        self.s_stdarray = np.zeros(self.arraylen)
        for evt in self.p_std.keys():
            for n in range(len(self.p_std[evt])):
                self.p_stdarray[index] = self.p_std[evt][n]
                self.s_stdarray[index] = self.s_std[evt][n]
                index += 1


    def histplot(self, num, label=None):
        binlist = []
        if num == 0:
            bfactor = 0.001
            badd = 0
        elif num == 1:
            bfactor = 0.003
            badd = 0
        else:
            bfactor = 0.006
            badd = 0.4
        for i in range(100):
            binlist.append(badd+bfactor*i)
        import matplotlib.pyplot as plt

        plt.hist(self.getTheta(num),bins = binlist)
        plt.xlabel('Values')
        plt.ylabel('Frequency')
        title_str = 'Quantile distance quotient distribution'
        if label:
            title_str += ' (' + label + ')'
        plt.title(title_str)
        plt.show()


    def getTheta(self,number):
        if number == 0:
            return self.theta0
        elif number == 1:
            return self.theta1
        elif number == 2:
            return self.theta2


    def getPDFDict(self, month, evt):
        self.pdfdict = PDFDictionary.from_quakeml(os.path.join(self.directory,month,evt))


    def getStatistics(self):
        self.p_mean = self.p_stdarray.mean()
        self.p_std_std = self.p_stdarray.std()
        self.p_median = np.median(self.p_stdarray)
        self.s_mean = self.s_stdarray.mean()
        self.s_std_std = self.s_stdarray.std()
        self.s_median = np.median(self.s_stdarray)


    def writeThetaToFile(self,array,out_dir,filename = None):
        fid = open(os.path.join(out_dir,filename), 'w')
        for val in array:
            fid.write(str(val)+'\n')
        fid.close()



#if __name__ == "__main__":