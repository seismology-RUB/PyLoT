#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import operator
import os

import matplotlib.pyplot as plt
import numpy as np
from obspy import read_events
from obspy.core import AttribDict
from pylot.core.io.phases import picksdict_from_picks
from pylot.core.util.pdf import ProbabilityDensityFunction
from pylot.core.util.utils import find_in_list
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
        self._pdfs = dict()
        names = self.iter_kwargs(kwargs)
        if len(names) > 2:
            raise ValueError('Comparison is only defined for two '
                             'arguments!')
        self._names = names
        self._compare = self.compare_picksets()

    def __nonzero__(self):
        if not len(self.names) == 2 or not self._pdfs:
            return False
        return True

    def iter_kwargs(self, kwargs):
        names = list()
        for name, fn in kwargs.items():
            if name == 'eventlist':
                names = self.init_by_eventlist(fn)
                break
            if isinstance(fn, PDFDictionary):
                self._pdfs[name] = fn
            elif isinstance(fn, dict) or isinstance(fn, AttribDict):
                self._pdfs[name] = PDFDictionary(fn)
            else:
                self._pdfs[name] = PDFDictionary.from_quakeml(fn)
            names.append(name)
        return names

    def init_by_eventlist(self, eventlist):
        # create one dictionary containing all picks for all events (therefore modify station key)
        global_picksdict = {}
        for event in eventlist:
            automanu = {'manu': event.pylot_picks,
                        'auto': event.pylot_autopicks}
            for method, picksdict in automanu.items():
                if not method in global_picksdict.keys():
                    global_picksdict[method] = {}
                for station, picks in picksdict.items():
                    new_picksdict = global_picksdict[method]
                    # new id combining event and station in one dictionary for all events
                    id = '{}_{}'.format(event.pylot_id, station)
                    new_picksdict[id] = picks
        for method, picksdict in global_picksdict.items():
            self._pdfs[method] = PDFDictionary(picksdict)
        names = list(global_picksdict.keys())
        return names

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

        pdf_a = self.get('auto').generate_pdf_data(type)
        pdf_b = self.get('manu').generate_pdf_data(type)

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
        method = operator.methodcaller(method_name)
        pdf_list = self.get_all(phase)
        rarray = list(map(method, pdf_list))
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
        return PDFDictionary(fn)

    def get_all(self, phase):
        rlist = list()
        for phases in self.pdf_data.values():
            try:
                rlist.append(phases[phase])
            except KeyError:
                continue
        return rlist

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
    """
    This object can be used to get various statistic values from probabillity density functions.
    Takes a path as argument.
    """

    def __init__(self, directory):
        """Initiates some values needed when dealing with pdfs later"""
        self._rootdir = directory
        self._evtlist = list()
        self._rphase = None
        self.make_fnlist()

    def make_fnlist(self, fn_pattern='*.xml'):
        """
        Takes a file pattern and searches for that recursively in the set path for the object.
        :param fn_pattern: A pattern that can identify all datafiles. Default Value = '*.xml'
        :type  fn_pattern: string
        :return: creates a list of events saved in the PDFstatistics object.
        """
        evtlist = list()
        for root, _, files in os.walk(self.root):
            for file in files:
                if file.endswith(fn_pattern[1:]):
                    evtlist.append(os.path.join(root, file))
        self._evtlist = evtlist

    def __iter__(self):
        for evt in self._evtlist:
            yield PDFDictionary.from_quakeml(evt)

    def __getitem__(self, item):
        evt = find_in_list(self._evtlist, item)
        if evt:
            return PDFDictionary.from_quakeml(evt)
        return None

    @property
    def root(self):
        return self._rootdir

    @root.setter
    def root(self, value):
        if os.path.exists(value):
            self._rootdir = value
        else:
            raise ValueError("path doesn't exist: %s" % value)

    @property
    def curphase(self):
        """
        return the current phase type of interest
        :return: current phase
        """
        return self._rphase

    @curphase.setter
    def curphase(self, type):
        """
        setter method for property curphase
        :param type: specify the phase type of interest
        :type  type: string ('p' or 's')
        :return: -
        """
        if type.upper() not in 'PS':
            raise ValueError("phase type must be either 'P' or 'S'!")
        else:
            self._rphase = type.upper()

    def get(self, property='std', value=None):
        """
        takes a property str and a probability value and returns all
        property's values for the current phase of interest
        :func:`self.curphase`

        :param property: property name (default: 'std')
        :type property: str
        :param value: probability value :math:\alpha
        :type  value: float
        :return: list containing all property's values
        """
        assert isinstance(self.curphase,
                          str), 'phase has to be set before being ' \
                                'able to iterate over items...'
        rlist = []
        method_options = dict(STD='standard_deviation',
                              Q='quantile',
                              QD='quantile_distance',
                              QDF='quantile_dist_frac')

        # create method caller for easy mapping
        if property.upper() == 'STD':
            method = operator.methodcaller(method_options[property.upper()])
        elif value is not None:
            try:
                method = operator.methodcaller(method_options[property.upper()],
                                               value)
            except KeyError:
                raise KeyError('unknwon property: {0}'.format(property.upper()))
        else:
            raise ValueError("for call to method {0} value has to be "
                             "defined but is 'None' ".format(method_options[
                                                                 property.upper()]))

        for pdf_dict in self:
            # create worklist
            wlist = pdf_dict.get_all(self.curphase)
            # map method calls to object in worklist
            rlist += map(method, wlist)

        return rlist

    def writeThetaToFile(self, array, out_dir):
        """
        Method to write array like data to file. Useful since acquiring can take
        serious amount of time when dealing with large databases.
        :param array: List of values.
        :type  array: list
        :param out_dir: Path to save file to including file name.
        :type  out_dir: str
        :return: Saves a file at given output directory.
        """
        fid = open(os.path.join(out_dir), 'w')
        for val in array:
            fid.write(str(val) + '\n')
        fid.close()


def main():
    root_dir = '/home/sebastianp/Codetesting/xmls/'
    Insheim = PDFstatistics(root_dir)
    Insheim.curphase = 'p'
    qdlist = Insheim.get('qdf', 0.2)
    print(qdlist)


if __name__ == "__main__":
    import cProfile

    pr = cProfile.Profile()
    pr.enable()
    main()
    pr.disable()
    # after your program ends
    pr.print_stats(sort="calls")
