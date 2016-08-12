#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hashlib
import numpy as np
from scipy.interpolate import splrep, splev
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

def fit_curve(x, y):

    return splev, splrep(x, y)

def getindexbounds(f, eta):
    mi = f.argmax()
    m = max(f)
    b = m * eta
    l = find_nearest(f[:mi], b)
    u = find_nearest(f[mi:], b) + mi
    return mi, l, u

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

def clims(lim1, lim2):
    """
    takes two pairs of limits and returns one pair of common limts
    :param lim1:
    :param lim2:
    :return:

    >>> clims([0, 4], [1, 3])
    [0, 4]
    >>> clims([1, 4], [0, 3])
    [0, 4]
    >>> clims([1, 3], [0, 4])
    [0, 4]
    >>> clims([0, 3], [1, 4])
    [0, 4]
    >>> clims([0, 3], [0, 4])
    [0, 4]
    >>> clims([1, 4], [0, 4])
    [0, 4]
    >>> clims([0, 4], [0, 4])
    [0, 4]
    >>> clims([0, 4], [1, 4])
    [0, 4]
    >>> clims([0, 4], [0, 3])
    [0, 4]
    """
    lim = [None, None]
    if lim1[0] < lim2[0]:
        lim[0] = lim1[0]
    else:
        lim[0] = lim2[0]
    if lim1[1] > lim2[1]:
        lim[1] = lim1[1]
    else:
        lim[1] = lim2[1]
    return lim


def demeanTrace(trace, window):
    """
    takes a trace object and returns the same trace object but with data
    demeaned within a certain time window
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
    :type combo_box: `~QComboBox`
    :param val: Name of a combo box to search for.
    :type val: basestring
    :return: index value of item with name val or 0
    """
    return combo_box.findText(val) if combo_box.findText(val) is not -1 else 0


def find_nearest(array, value):
    '''
    function find_nearest takes an array and a value and returns the
    index of the nearest value found in the array
    :param array: array containing values
    :type array: `~numpy.ndarray`
    :param value: number searched for
    :return: index of the array item being nearest to the value

    >>> a = np.array([ 1.80339578, -0.72546654,  0.95769195, -0.98320759, 0.85922623])
    >>> find_nearest(a, 1.3)
    2
    >>> find_nearest(a, 0)
    1
    >>> find_nearest(a, 2)
    0
    >>> find_nearest(a, -1)
    3
    >>> a = np.array([ 1.1, -0.7,  0.9, -0.9, 0.8])
    >>> find_nearest(a, 0.849)
    4
    '''
    return (np.abs(array - value)).argmin()


def fnConstructor(s):
    '''
    takes a string and returns a valid filename (especially on windows machines)
    :param s: desired filename
    :type s: str
    :return: valid filename
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


def four_digits(year):
    """
    takes a two digit year integer and returns the correct four digit equivalent
    from the last 100 years
    :param year: two digit year
    :type year: int
    :return: four digit year correspondant

    >>> four_digits(20)
    1920
    >>> four_digits(16)
    2016
    >>> four_digits(00)
    2000
    """
    if year + 2000 <= UTCDateTime.utcnow().year:
        year += 2000
    else:
        year += 1900
    return year


def getGlobalTimes(stream):
    '''
    takes a stream object and returns the latest end and the earliest start
    time of all contained trace objects
    :param stream: seismological data stream
    :type stream: `~obspy.core.stream.Stream`
    :return: minimum start time and maximum end time
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
    takes a time object and returns the corresponding SHA1 hash of the
    formatted date string
    :param time: time object for which a hash should be calculated
    :type time: :class: `~obspy.core.utcdatetime.UTCDateTime` object
    :return: str
    '''
    hg = hashlib.sha1()
    hg.update(time.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return hg.hexdigest()


def getLogin():
    '''
    returns the actual user's login ID
    :return: login ID
    '''
    return pwd.getpwuid(os.getuid())[0]


def getOwner(fn):
    '''
    takes a filename and return the login ID of the actual owner of the file
    :param fn: filename of the file tested
    :type fn: str
    :return: login ID of the file's owner
    '''
    return pwd.getpwuid(os.stat(fn).st_uid).pw_name


def getPatternLine(fn, pattern):
    """
    takes a file name and a pattern string to search for in the file and
    returns the first line which contains the pattern string otherwise 'None'
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
    takes an iterable and returns 'True' if the items are in order otherwise
    'False'
    :param iterable: an iterable object
    :type iterable:
    :return: Boolean

    >>> isSorted(1)
    Traceback (most recent call last):
    ...
    AssertionError: object is not iterable; object: 1
    >>> isSorted([1,2,3,4])
    True
    >>> isSorted('abcd')
    True
    >>> isSorted('bcad')
    False
    >>> isSorted([2,3,1,4])
    False
    '''
    assert isIterable(iterable), 'object is not iterable; object: {' \
                                 '0}'.format(iterable)
    if type(iterable) is str:
        iterable = [s for s in iterable]
    return sorted(iterable) == iterable


def isIterable(obj):
    """
    takes a python object and returns 'True' is the object is iterable and
    'False' otherwise
    :param obj: a python object
    :return: True of False
    """
    try:
        iterator = iter(obj)
    except TypeError as te:
        return False
    return True

def prepTimeAxis(stime, trace):
    '''
    takes a starttime and a trace object and returns a valid time axis for
    plotting
    :param stime: start time of the actual seismogram as UTCDateTime
    :param trace: seismic trace object
    :return: valid numpy array with time stamps for plotting
    '''
    nsamp = trace.stats.npts
    srate = trace.stats.sampling_rate
    tincr = trace.stats.delta
    etime = stime + nsamp / srate
    time_ax = np.arange(stime, etime, tincr)
    if len(time_ax) < nsamp:
        print('elongate time axes by one datum')
        time_ax = np.arange(stime, etime + tincr, tincr)
    elif len(time_ax) > nsamp:
        print('shorten time axes by one datum')
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
