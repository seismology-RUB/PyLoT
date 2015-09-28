#!/usr/bin/env python
# -*- coding: utf-8 -*-


class AutoPickParameter(object):
    '''
    AutoPickParameters is a parameter type object capable to read and/or write
    parameter ASCII.

    :param fn str: Filename of the input file

    Parameters are given for example as follows:
    ==========  ==========  =======================================
    Name        Value       Comment
    ==========  ==========  =======================================
    phl         S           # phaselabel
    ff1         0.1         # freqmin
    ff2         0.5         # freqmax
    tdet        6.875       # det-window_(s)_for_ar
    tpred       2.5         # pred-window_(s)_for_ar
    order       4           # order_of_ar
    fnoise      0           # noise_level_for_ar
    suppp       7           # envelopecoeff
    tolt        300         # (s)time around arrival time
    f1tpwt      4           # propfact_minfreq_secondtaper
    pickwindow  9           # length_of_pick_window
    w1          1           # length_of_smoothing_window
    w2          0.37        # cf(i-1)*(1+peps)_for_local_min
    w3          0.25        # cf(i-1)*(1+peps)_for_local_min
    tslope      0.8;2       # slope_det_window_loc_glob
    aerr        30;60       # adjusted_error_slope_fitting_loc_glob
    tsn         20;5;20;10  # length_signal_window_S/N
    proPh       Sn          # nextprominentphase
    ==========  ==========  =======================================
    '''

    def __init__(self, fnin=None, fnout=None, **kwargs):
        '''
        Initialize parameter object:

        read content of an ASCII file an form a type consistent dictionary
        contain all parameters.
        '''

        self.__filename = fnin
        parFileCont = {}
        # read from parsed arguments alternatively
        for key, val in kwargs.items():
            parFileCont[key] = val

        if self.__filename is not None:
                inputFile = open(self.__filename, 'r')
        else:
            return
        try:
            lines = inputFile.readlines()
            for line in lines:
                parspl = line.split('\t')[:2]
                parFileCont[parspl[0].strip()] = parspl[1]
        except Exception as e:
            self._printParameterError(e)
            inputFile.seek(0)
            lines = inputFile.readlines()
            for line in lines:
                if not line.startswith(('#', '%', '\n', ' ')):
                    parspl = line.split('#')[:2]
                    parFileCont[parspl[1].strip()] = parspl[0].strip()
        for key, value in parFileCont.items():
            try:
                val = int(value)
            except:
                try:
                    val = float(value)
                except:
                    if len(value.split(' ')) > 1:
                        vallist = value.strip().split(' ')
                        val = []
                        for val0 in vallist:
                            val0 = float(val0)
                            val.append(val0)
                    else:
                        val = str(value.strip())
            parFileCont[key] = val
        self.__parameter = parFileCont

        if fnout:
            self.export2File(fnout)

    # Human-readable string representation of the object
    def __str__(self):
        string = ''
        string += 'Automated picking parameter:\n\n'
        if self.__parameter:
            for key, value in self.iteritems():
                string += '%s:\t\t%s\n' % (key, value)
        else:
            string += 'Empty parameter dictionary.'
        return string

    # String representation of the object
    def __repr__(self):
        return "AutoPickParameter('%s')" % self.__filename

    # Boolean test
    def __nonzero__(self):
        return self.__parameter

    def __getitem__(self, key):
        return self.__parameter[key]

    def __setitem__(self, key, value):
        self.__parameter[key] = value

    def __delitem__(self, key):
        del self.__parameter[key]

    def __iter__(self):
        return iter(self.__parameter)

    def __len__(self):
        return len(self.__parameter.keys())

    def iteritems(self):
        for key, value in self.__parameter.items():
            yield key, value

    def hasParam(self, parameter):
        if self.__parameter.has_key(parameter):
            return True
        return False

    def getParam(self, *args):
        try:
            for param in args:
                try:
                    return self.__getitem__(param)
                except KeyError as e:
                    self._printParameterError(e)
        except TypeError:
            try:
                return self.__getitem__(args)
            except KeyError as e:
                self._printParameterError(e)

    def setParam(self, **kwargs):
        for param, value in kwargs.items():
            self.__setitem__(param, value)
        print(self)

    @staticmethod
    def _printParameterError(errmsg):
        print('ParameterError:\n non-existent parameter %s' % errmsg)

    def export2File(self, fnout):
        fid_out = open(fnout, 'w')
        lines = []
        for key, value in self.iteritems():
            lines.append('{key}\t{value}'.format(key=key, value=value))
        fid_out.writelines(lines)


class FilterOptions(object):
    '''
    FilterOptions is a parameter object type providing Butterworth filter
    option parameter for PyLoT. Its easy to access properties helps to manage
    file based as well as parameter manipulation within the GUI.

    :type filtertype: str, optional
    :param filtertype: String containing the desired filtertype For information
    about the supported filter types see _`Supported Filter` section .

    :type freq: list, optional
    :param freq: list of float(s) describing the cutoff limits of the filter

    :type order: int, optional
    :param order: Integer value describing the order of the desired Butterworth
    filter.

    .. rubric:: _`Supported Filter`

        ``'bandpass'``
            Butterworth-Bandpass

        ``'bandstop'``
            Butterworth-Bandstop

        ``'lowpass'``
            Butterworth-Lowpass

        ``'highpass'``
            Butterworth-Highpass
    '''
    def __init__(self, filtertype='bandpass', freq=[2., 5.], order=3,
                 **kwargs):
        self._order = order
        self._filtertype = filtertype
        self._freq = freq

    def __str__(self):
        hrs = '''\n\tFilter parameter:\n
        Type:\t\t{ftype}\n
        Frequencies:\t{freq}\n
        Order:\t\t{order}\n
        '''.format(ftype=self.getFilterType(),
                   freq=self.getFreq(),
                   order=self.getOrder())
        return hrs

    def parseFilterOptions(self):
        if self.getFilterType():
            robject = {'type':self.getFilterType()}
            robject['corners'] = self.getOrder()
            if len(self.getFreq()) > 1:
                robject['freqmin'] = self.getFreq()[0]
                robject['freqmax'] = self.getFreq()[1]
            else:
                robject['freq'] = self.getFreq() if type(self.getFreq()) is \
                                                    float else self.getFreq()[0]
            return robject
        return None

    def getFreq(self):
        return self.__getattribute__('_freq')

    def setFreq(self, freq):
        self.__setattr__('_freq', freq)

    def getOrder(self):
        return self.__getattribute__('_order')

    def setOrder(self, order):
        self.__setattr__('_order', order)

    def getFilterType(self):
        return self.__getattribute__('_filtertype')

    def setFilterType(self, filtertype):
        self.__setattr__('_filtertype', filtertype)
