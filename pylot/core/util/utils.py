#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hashlib
import numpy as np
import os
import pwd
import re
import subprocess
from obspy.core import UTCDateTime


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)


def worker(func, input, cores='max', async=False):
    import multiprocessing

    if cores == 'max':
        cores = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(cores)

    if async == True:
        result = pool.map_async(func, input)
    else:
        result = pool.map(func, input)
    pool.close()
    return result


def demeanTrace(trace, window):
    """
    returns the DATA where each trace is demean by the average value within
    WINDOW
    :param trace: waveform trace object
    :type trace: `~obspy.core.stream.Trace`
    :param window:
    :type window: tuple
    :return: trace
    :rtype: `~obspy.core.stream.Trace`
    """
    trace.data -= trace.data[window].mean()
    return trace


def findComboBoxIndex(combo_box, val):
    """
    Function findComboBoxIndex takes a QComboBox object and a string and
    returns either 0 or the index throughout all QComboBox items.
    :param combo_box: Combo box object.
    :type combo_box: QComboBox
    :param val: Name of a combo box to search for.
    :type val:
    :return: index value of item with name val or 0
    """
    return combo_box.findText(val) if combo_box.findText(val) is not -1 else 0


def find_nearest(array, value):
    '''
    Function find_nearest takes an array and a value and returns the
    index of the nearest value found in the array.
    :param array:
    :param value:
    :return:
    '''
    return (np.abs(array - value)).argmin()


def fnConstructor(s):
    '''

    :param s:
    :type s:
    :return:
    '''
    if type(s) is str:
        s = s.split(':')[-1]
    else:
        s = getHash(UTCDateTime())

    badchars = re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
    badsuffix = re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

    fn = badchars.sub('_', s)

    if badsuffix.match(fn):
        fn = '_' + fn
    return fn


def getGlobalTimes(stream):
    '''

    :param stream:
    :type stream
    :return:
    '''
    min_start = UTCDateTime()
    max_end = None
    for trace in stream:
        if trace.stats.starttime < min_start:
            min_start = trace.stats.starttime
        if max_end is None or trace.stats.endtime > max_end:
            max_end = trace.stats.endtime
    return min_start, max_end


def getHash(time):
    '''
    :param time: time object for which a hash should be calculated
    :type time: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :return: str
    '''
    hg = hashlib.sha1()
    hg.update(time.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return hg.hexdigest()


def getLogin():
    '''

    :return:
    '''
    return pwd.getpwuid(os.getuid())[0]


def getOwner(fn):
    '''

    :param fn:
    :type fn:
    :return:
    '''
    return pwd.getpwuid(os.stat(fn).st_uid).pw_name


def getPatternLine(fn, pattern):
    """
    Takes a file name and a pattern string to search for in the file and
    returns the first line which contains the pattern string otherwise None.

    :param fn: file name
    :type fn: str
    :param pattern: pattern string to search for
    :type pattern: str
    :return: the complete line containing pattern or None

    >>> getPatternLine('utils.py', 'python')
    '#!/usr/bin/env python\\n'
    >>> print(getPatternLine('version.py', 'palindrome'))
    None
    """
    fobj = open(fn, 'r')
    for line in fobj.readlines():
        if pattern in line:
            fobj.close()
            return line

    return None


def isSorted(iterable):
    '''

    :param iterable:
    :type iterable:
    :return:
    '''
    return sorted(iterable) == iterable


def prepTimeAxis(stime, trace):
    '''

    :param stime:
    :type stime:
    :param trace:
     :type trace:
    :return:
    '''
    nsamp = trace.stats.npts
    srate = trace.stats.sampling_rate
    tincr = trace.stats.delta
    etime = stime + nsamp / srate
    time_ax = np.arange(stime, etime, tincr)
    if len(time_ax) < nsamp:
        print 'elongate time axes by one datum'
        time_ax = np.arange(stime, etime + tincr, tincr)
    elif len(time_ax) > nsamp:
        print 'shorten time axes by one datum'
        time_ax = np.arange(stime, etime - tincr, tincr)
    if len(time_ax) != nsamp:
        raise ValueError('{0} samples of data \n '
                         '{1} length of time vector \n'
                         'delta: {2}'.format(nsamp, len(time_ax), tincr))
    return time_ax


def scaleWFData(data, factor=None, components='all'):
    """
    produce scaled waveforms from given waveform data and a scaling factor,
    waveform may be selected by their components name
    :param data: waveform data to be scaled
    :type data: `~obspy.core.stream.Stream` object
    :param factor: scaling factor
    :type factor: float
    :param components: components labels for the traces in data to be scaled by
     the scaling factor (optional, default: 'all')
    :type components: tuple
    :return:  scaled waveform data
    :rtype: `~obspy.core.stream.Stream` object
    """
    if components is not 'all':
        for comp in components:
            if factor is None:
                max_val = np.max(np.abs(data.select(component=comp)[0].data))
                data.select(component=comp)[0].data /= 2 * max_val
            else:
                data.select(component=comp)[0].data /= 2 * factor
    else:
        for tr in data:
            if factor is None:
                max_val = float(np.max(np.abs(tr.data)))
                tr.data /= 2 * max_val
            else:
                tr.data /= 2 * factor

    return data


def runProgram(cmd, parameter=None):
    """
    run an external program specified by cmd with parameters input returning the
    stdout output

    :param cmd: name of the command to run
    :type cmd: str
    :param parameter: filename of parameter file  or parameter string
    :type parameter: str
    :return: stdout output
    :rtype: str
    """

    if parameter:
        cmd.strip()
        cmd += ' %s 2>&1' % parameter

    output = subprocess.check_output('{} | tee /dev/stderr'.format(cmd),
                                     shell=True)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
