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
    """
    This object can be used to get various statistic values from probabillity density functions.
    Takes a path as argument.
    """


    def __init__(self, directory):
        """Initiates some values needed when dealing with pdfs later"""
        self.directory = directory
        self.evtlist = list()
        self.return_phase = None


    def readTheta(self, arname, dir, fnpattern):
        """
        Loads an array from file into object instance.
        :param arname: Name of Array beeing created.
        :type  arname: string
        :param dir: Directory where file is to be found.
        :type  dir: string
        :param fnpattern: file name pattern for reading multiple files into one array.
        :type  fnpattern: string
        :return: a list with all args* from the files.
        """
        exec('self.' + arname +' = []')
        filelist = glob.glob1(dir, fnpattern)
        for file in filelist:
            fid = open(os.path.join(dir,file), 'r')
            list = []
            for line in fid.readlines():
                list.append(eval(line))
            exec('self.' + arname + ' += list')
            fid.close()


    def makeFileList(self, fn_pattern='*.xml'):
        """
        Takes a file pattern and searches for that recursively in the set path for the object.
        :param fn_pattern: A pattern that can identify all datafiles. Default Value = '*.xml'
        :type  fn_pattern: string
        :return: creates a list of events saved in the PDFstatistics object.
        """
        evtlist = glob.glob1((os.path.join(self.directory)), fn_pattern)
        if not evtlist:
             for root, _, files in os.walk(self.directory):
                 for file in files:
                     if file.endswith(fn_pattern[1:]):
                         evtlist.append(os.path.join(root, file))
        self.evtlist = evtlist


    def __iter__(self):
        """Iterating over the PDFstatistics object yields every single pdf from the list of events"""
        assert isinstance(self.return_phase, str), 'phase has to be set before being able to iterate over items...'
        for evt in self.evtlist:
            self.getPDFDict(self.directory, evt)
            for station, pdfs in self.pdfdict.pdf_data.items():
                try:
                    yield pdfs[self.return_phase]
                except KeyError:
                    continue


    def set_return_phase(self, type):
        """
        Sets the phase typ of event data that is returned on iteration over the object.
        :param type: can be either p (p-phase) or s (s-phase).
        :type  type: string
        :return: -
        """
        if type.upper() not in 'PS':
            raise ValueError("phase type must be either 'P' or 'S'!")
        else:
            self.return_phase = type.upper()


    def getQD(self,value):
        """
        Takes a probability value and and returns the distance
        between two complementary quantiles.
        For example: getQD(0.3) yields Quantile(1-0.3) - Quantile(0.3)
        :param value: 0 < value < 0.5
        :type  value: float
        :return: returns a list of all quantile distances for all pdfs in
                 the list of events.
        """
        QDlist = []
        for pdf in self:
            QD = pdf.quantile_distance(value)
            QDlist.append(QD)
        return QDlist


    def getQDQ(self,value):
        """
        Takes a probability value and and returns the fraction of
        two quantile distances.
        For example:
        getQDQ(x) = getQD(0.5-x)/getQD(x)
        (Quantile(1-0.5-x) - Quantile(x)) / (Quantile(1-x) - Quantile(x))
        :param value: 0 < value < 0.25
        :return: returns a list of all quantile fractions for all pdfs in
                 the list of events.
        """
        QDQlist = []
        for pdf in self:
            QDQ = pdf.qtile_dist_quot(value)
            QDQlist.append(QDQ)
        return QDQlist


    def getSTD(self):
        """
        Iterates over PDFstatistics object and returns the standard
        deviation of all pdfs in the list of events.
        :return: saves an instance of self.p_stdarray or
                 self.s_stdarray, depending on set phase.
        """
        std = []
        for pdf in self:
            try:
                std.append(pdf.standard_deviation())
            except KeyError:
                continue
        std = np.array(std)
        self.set_stdarray(std)


    def set_stdarray(self, array):
        """
        Helper function for self.getSTD(). This function
        should not be called directly.
        """
        if self.return_phase == 'P':
            self.p_stdarray = array
        elif self.return_phase == 'S':
            self.s_stdarray = array
        else:
            raise ValueError('phase type not set properly...\n'
                             'Actual phase type: {0}'.format(self.return_phase))


    def getBinList(self,l_boundary,u_boundary,nbins = 100):
        """
        Helper function for self.histplot(). Takes in two boundaries and
        a number of bins and creates a list of bins which can be passed
        to self.histplot().
        :param l_boundary: Any number.
        :type  l_boundary: float
        :param u_boundary: Any number that is greater than l_boundary.
        :type  u_boundary: float
        :param nbins: Any positive integer.
        :type  nbins: int
        :return: A list of equidistant bins.
        """
        if u_boundary <= l_boundary:
            raise ValueError('Upper boundary must be greather than lower!')
        elif nbins <= 0:
            raise ValueError('Number of bins is not valid.')
        binlist = []
        for i in range(nbins):
            binlist.append(l_boundary + i*(u_boundary-l_boundary)/nbins)
        return  binlist


    def histplot(self, array, binlist, xlab = 'Values',
                 ylab = 'Frequency', title = None, fnout = None):
        """
        Method to quickly show some distribution of data. Takes array like data,
        and a list of bins. Editing detail and inserting a legend is not possible.
        :param array: List of values.
        :type  array: Array like
        :param binlist: List of bins.
        :type  binlist: list
        :param xlab: A label for the x-axes.
        :type  xlab: str
        :param ylab: A label for the y-axes.
        :type  ylab: str
        :param title: A title for the Plot.
        :type  title: str
        :param fnout: A path to save the plot instead of showing.
                      Has to contain filename and type. Like: 'path/to/file.png'
        :type  fnout. str
        :return: -
        """
        import matplotlib.pyplot as plt
        plt.hist(array,bins = binlist)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        if title:
            plt.title(title)
        if fnout:
            plt.savefig(fnout)
        else:
            plt.show()


    def getPDFDict(self, month, evt):
        """
        Helper function for __iter__(). Should not be called directly.
        """
        self.pdfdict = PDFDictionary.from_quakeml(os.path.join(self.directory,month,evt))


    def getStatistics(self):
        """
        On call function will get mean, median and standard deviation values
        from self.p_stdarray and self.s_stdarray. Both must be
        instances before calling this function.
        :return: Creates instances of self.p_mean, self.p_std_std and self.p_median
                 for both phases (simultaneously) for the PDFstatistics object.
        """
        if not self.p_stdarray or not self.s_stdarray:
            raise NotImplementedError('Arrays are not properly set yet!')
        elif type(self.p_stdarray) != type(np.zeros(1)) or type(self.s_stdarray) != type(np.zeros(1)):
            raise TypeError('Array is not a proper numpy array.')

        self.p_mean = self.p_stdarray.mean()
        self.p_std_std = self.p_stdarray.std()
        self.p_median = np.median(self.p_stdarray)
        self.s_mean = self.s_stdarray.mean()
        self.s_std_std = self.s_stdarray.std()
        self.s_median = np.median(self.s_stdarray)


    def writeThetaToFile(self,array,out_dir):
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
            fid.write(str(val)+'\n')
        fid.close()


def main():
    root_dir ='/home/sebastianp/Codetesting/xmls/'
    Insheim = PDFstatistics(root_dir)
    Insheim.makeFileList()
    Insheim.set_return_phase('p')
    Insheim.getSTD()
    qdlist = Insheim.getQDQ(0.3)
    binlist = Insheim.getBinList(0.,3.)
    print qdlist


if __name__ == "__main__":
    main()