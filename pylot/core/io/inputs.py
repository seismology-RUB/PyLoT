#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import os

from pylot.core.io import default_parameters
from pylot.core.util.errors import ParameterError


class PylotParameter(object):
    '''
    PylotParameter is a parameter type object capable to read and/or write
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

    def __init__(self, fnin=None, fnout=None, verbosity=0, **kwargs):
        '''
        Initialize parameter object:

        io content of an ASCII file an form a type consistent dictionary
        contain all parameters.
        '''

        self.__init_default_paras()
        self.__init_subsettings()
        self.__filename = fnin
        self.__parameter = {}
        self._verbosity = verbosity
        self._parFileCont = {}
        # io from parsed arguments alternatively
        for key, val in kwargs.items():
            self._parFileCont[key] = val
        self.from_file()
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

    # Set default values of parameter names
    def __init_default_paras(self):
        """set default values of parameter names"""
        parameters = default_parameters.defaults
        self.__defaults = parameters

    def __init_subsettings(self):
        self._settings_main = default_parameters.settings_main
        self._settings_special_pick = default_parameters.settings_special_pick

    # String representation of the object
    def __repr__(self):
        return "PylotParameter('%s')" % self.__filename

    # Boolean test
    def __nonzero__(self):
        return bool(self.__parameter)

    def __getitem__(self, key):
        if key in self.__parameter:
            return self.__parameter[key]
        else:
            logging.warning(f'{key} not found in PylotParameter')

    def __setitem__(self, key, value):
        try:
            value = self.check_range(value, self.__defaults[key]['max'], self.__defaults[key]['min'])
        except KeyError:
            # no min/max values in defaults
            pass
        self.__parameter[key] = value

    def __delitem__(self, key):
        del self.__parameter[key]

    def __iter__(self):
        return iter(self.__parameter)

    def __len__(self):
        return len(self.__parameter.keys())

    def iteritems(self):
        """
        Iterate over parameters
        :return: key, value tupel
        :rtype:
        """
        for key, value in self.__parameter.items():
            yield key, value

    def hasParam(self, parameter):
        """
        Check if parameter is in keys
        :param parameter: parameter to look for in keys
        :type parameter:
        :return:
        :rtype: bool
        """
        return parameter in self.__parameter.keys()

    def get(self, *args):
        """
        Get first available parameter in args
        :param args:
        :type args:
        :return:
        :rtype:
        """
        try:
            for param in args:
                try:
                    return self.__getitem__(param)
                except KeyError as e:
                    self._printParameterError(e)
                    raise ParameterError(e)
        except TypeError:
            try:
                return self.__getitem__(args)
            except KeyError as e:
                self._printParameterError(e)
                raise ParameterError(e)

    def get_defaults(self):
        """
        get default parameters
        :return:
        :rtype: dict
        """
        return self.__defaults

    def get_main_para_names(self):
        """
        Get main parameter names
        :return: list of keys available in parameters
        :rtype:
        """
        return self._settings_main

    def get_special_para_names(self):
        """
        Get pick parameter names
        :return: list of keys available in parameters
        :rtype:
        """
        return self._settings_special_pick

    def get_all_para_names(self):
        """
        Get all parameter names
        :return:
        :rtype: list
        """
        all_names = []
        all_names += self.get_main_para_names()['dirs']
        all_names += self.get_main_para_names()['nlloc']
        all_names += self.get_main_para_names()['smoment']
        all_names += self.get_main_para_names()['localmag']
        all_names += self.get_main_para_names()['pick']
        all_names += self.get_main_para_names()['filter']
        all_names += self.get_special_para_names()['z']
        all_names += self.get_special_para_names()['h']
        all_names += self.get_special_para_names()['fm']
        all_names += self.get_special_para_names()['quality']
        return all_names

    def reinit_default_parameters(self):
        self.__init_default_paras()

    @staticmethod
    def check_range(value, max_value, min_value):
        """
        Check if value is within the min/max values defined in default_parameters. Works for tuple and scalar values.
        :param value: Value to be checked against min/max range
        :param max_value: Maximum allowed value, tuple or scalar
        :param min_value: Minimum allowed value, tuple or scalar
        :return: value tuple/scalar clamped to the valid range

        >>> checkRange(-5, 10, 0)
         0
         >>> checkRange((-5., 100.), (10., 10.), (0., 0.))
         (0.0, 10.0)
        """
        try:
            # Try handling tuples by comparing their elements
            comparisons = [(a > b) for a, b in zip(value, max_value)]
            if True in comparisons:
                value = tuple(max_value[i] if comp else value[i] for i, comp in enumerate(comparisons))
            comparisons = [(a < b) for a, b in zip(value, min_value)]
            if True in comparisons:
                value = tuple(min_value[i] if comp else value[i] for i, comp in enumerate(comparisons))
        except TypeError:
            value = max(min_value, min(max_value, value))
        return value

    def checkValue(self, param, value):
        """
        Check type of value against expected type of param.
        Print warning message if type check fails
        :param param:
        :type param:
        :param value:
        :type value:
        :return:
        :rtype:
        """
        is_type = type(value)
        expect_type = self.get_defaults()[param]['type']
        if not is_type == expect_type and not is_type == tuple:
            message = 'Type check failed for param: {}, is type: {}, expected type:{}'
            message = message.format(param, is_type, expect_type)
            print(Warning(message))

    def setParamKV(self, param, value):
        """
        set parameter param to value
        :param param:
        :type param:
        :param value:
        :type value:
        :return:
        :rtype: None
        """
        self.__setitem__(param, value)

    def setParam(self, **kwargs):
        """
        Set multiple parameters
        :param kwargs:
        :type kwargs:
        :return:
        :rtype: None
        """
        for key in kwargs:
            self.__setitem__(key, kwargs[key])

    @staticmethod
    def _printParameterError(errmsg):
        print('ParameterError:\n non-existent parameter %s' % errmsg)

    def reset_defaults(self):
        """
        Reset current parameters to default parameters
        :return:
        :rtype: None
        """
        defaults = self.get_defaults()
        for param_name, param in defaults.items():
            self.setParamKV(param_name, param['value'])

    def from_file(self, fnin=None):
        """
        read parameters from file and set values to read values
        :param fnin: filename
        :type fnin:
        :return:
        :rtype: None
        """
        if not fnin:
            if self.__filename is not None:
                fnin = self.__filename
            else:
                return

        if isinstance(fnin, (list, tuple)):
            fnin = fnin[0]

        inputFile = open(fnin, 'r')
        try:
            lines = inputFile.readlines()
            for line in lines:
                parspl = line.split('\t')[:2]
                self._parFileCont[parspl[0].strip()] = parspl[1]
        except IndexError as e:
            if self._verbosity > 0:
                self._printParameterError(e)
            inputFile.seek(0)
            lines = inputFile.readlines()
            for line in lines:
                if not line.startswith(('#', '%', '\n', ' ')):
                    parspl = line.split('#')[:2]
                    self._parFileCont[parspl[1].strip()] = parspl[0].strip()
        for key, value in self._parFileCont.items():
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
                            try:
                                val0 = float(val0)
                            except:
                                pass
                            val.append(val0)
                    else:
                        val = str(value.strip())
            self._parFileCont[key] = val
        self.__parameter = self._parFileCont

    def export2File(self, fnout):
        """
        Export parameters to file
        :param fnout: Filename of export file
        :type fnout: str
        :return:
        :rtype:
        """
        fid_out = open(fnout, 'w')
        lines = []
        # for key, value in self.iteritems():
        #     lines.append('{key}\t{value}\n'.format(key=key, value=value))
        # fid_out.writelines(lines)
        header = ('%This is a parameter input file for PyLoT/autoPyLoT.\n' +
                  '%All main and special settings regarding data handling\n' +
                  '%and picking are to be set here!\n' +
                  '%Parameters are optimized for %{} data sets!\n'.format(self.get_main_para_names()['pick'][0]))
        separator = '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'

        fid_out.write(header)
        self.write_section(fid_out, self.get_main_para_names()['dirs'],
                           'main settings', separator)
        self.write_section(fid_out, self.get_main_para_names()['nlloc'],
                           'NLLoc settings', separator)
        self.write_section(fid_out, self.get_main_para_names()['smoment'],
                           'parameters for seismic moment estimation', separator)
        self.write_section(fid_out, self.get_main_para_names()['localmag'],
                           'settings local magnitude', separator)
        self.write_section(fid_out, self.get_main_para_names()['filter'],
                           'filter settings', separator)
        self.write_section(fid_out, self.get_main_para_names()['pick'],
                           'common settings picker', separator)
        fid_out.write(('#special settings for calculating CF#\n' +
                       '%!!Edit the following only if you know what you are doing!!%\n'))
        self.write_section(fid_out, self.get_special_para_names()['z'],
                           'Z-component', None)
        self.write_section(fid_out, self.get_special_para_names()['h'],
                           'H-components', None)
        self.write_section(fid_out, self.get_special_para_names()['fm'],
                           'first-motion picker', None)
        self.write_section(fid_out, self.get_special_para_names()['quality'],
                           'quality assessment', None)

    def write_section(self, fid, names, title, separator):
        """
        write a section of parameters to file
        :param fid: File object to write to
        :type fid:
        :param names: which parameter names to write to file
        :type names:
        :param title: title of section
        :type title: str
        :param separator: section separator, written at start of section
        :type separator: str
        :return:
        :rtype:
        """
        if separator:
            fid.write(separator)
        fid.write('#{}#\n'.format(title))
        l_val = 50
        l_name = 15
        l_ttip = 100
        for name in names:
            value = self[name]
            if type(value) == list or type(value) == tuple:
                value_tmp = ''
                for vl in value:
                    value_tmp += '{} '.format(vl)
                value = value_tmp
            tooltip = self.get_defaults()[name]['tooltip']
            if not len(str(value)) > l_val:
                value = '{:<{}} '.format(str(value), l_val)
            else:
                value = '{} '.format(str(value))
            name += '#'
            if not len(name) > l_name:
                name = '#{:<{}} '.format(name, l_name)
            else:
                name = '#{} '.format(name)
            if not len(tooltip) > l_ttip:
                ttip = '%{:<{}}\n'.format(tooltip, l_ttip)
            else:
                ttip = '%{}\n'.format(tooltip)
            line = value + name + ttip
            fid.write(line)

    @staticmethod
    def check_deprecated_parameters(parameters):
        if parameters.hasParam('database') and parameters.hasParam('rootpath'):
            parameters['datapath'] = os.path.join(parameters['rootpath'], parameters['datapath'],
                                                  parameters['database'])
            logging.warning(
                f'Parameters database and rootpath are deprecated. '
                f'Tried to merge them to now path: {parameters["datapath"]}.'
            )

        remove_keys = []
        for key in parameters:
            if not key in default_parameters.defaults.keys():
                remove_keys.append(key)
                logging.warning(f'Removing deprecated parameter: {key}')

        for key in remove_keys:
            del parameters[key]

        parameters._settings_main = default_parameters.settings_main
        parameters._settings_special_pick = default_parameters.settings_special_pick


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

    def __nonzero__(self):
        return bool(self.getFilterType())

    def parseFilterOptions(self):
        if self:
            robject = {'type': self.getFilterType(),
                       'corners': self.getOrder(),
                       'zerophase': False}
            if not self.getFilterType() in ['highpass', 'lowpass']:
                robject['freqmin'] = self.getFreq()[0]
                robject['freqmax'] = self.getFreq()[1]
            elif self.getFilterType() == 'highpass':
                robject['freq'] = self.getFreq()[0]
            elif self.getFilterType() == 'lowpass':
                robject['freq'] = self.getFreq()[1]
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
