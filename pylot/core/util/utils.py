#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hashlib
import numpy as np
import os
import platform
import re
import subprocess
import warnings
from obspy import UTCDateTime, read
from obspy.core import AttribDict
from obspy.signal.rotate import rotate2zne
from scipy.interpolate import splrep, splev

from pylot.core.io.inputs import PylotParameter, FilterOptions
from pylot.core.util.obspyDMT_interface import check_obspydmt_eventfolder
from pylot.styles import style_settings


def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)


def getAutoFilteroptions(phase, parameter):
    filtername = {'P': 'bpz2',
                  'S': 'bph2'}
    if not phase in filtername.keys():
        print('autoPickParameter: No filter options for phase {}.'.format(phase))
        return
    freqmin, freqmax = parameter.get(filtername[phase])
    filteroptions = FilterOptions(type='bandpass', freq=[freqmin, freqmax], order=4)  # order=4 default from obspy
    return filteroptions


def readDefaultFilterInformation(fname):
    """
    Read default filter information from pylot.in file
    :param fname: path to pylot.in file
    :type fname: str
    :return: dictionary containing the defailt filter information
    :rtype: dict
    """
    pparam = PylotParameter(fname)
    return readFilterInformation(pparam)


def readFilterInformation(pylot_parameter):
    """
    Read filter information from PylotParameter object into a dictionary
    :param pylot_parameter: PylotParameter object
    :type pylot_parameter: `~pylot.pylot.core.io.inputs.PylotParameter`
    :return: dictionary containing the filter information
    :rtype: dict
    """
    p_filter = {'filtertype': pylot_parameter['filter_type'][0],
                'freq': [pylot_parameter['minfreq'][0], pylot_parameter['maxfreq'][0]],
                'order': int(pylot_parameter['filter_order'][0])}
    s_filter = {'filtertype': pylot_parameter['filter_type'][1],
                'freq': [pylot_parameter['minfreq'][1], pylot_parameter['maxfreq'][1]],
                'order': int(pylot_parameter['filter_order'][1])}
    filter_information = {'P': p_filter,
                          'S': s_filter}
    return filter_information


def fit_curve(x, y):
    """

    :param x: data points defining a curve y = f(x)
    :type x: array_like
    :param y: data points defining a curve y = f(x)
    :type y: array_like
    :return: tuple containing a function to evaluate a B-spline and
    (the vector of knots, B-spline coefficients, degree of the spline)
    :rtype: (func, (t, c, k))
    """
    return splev, splrep(x, y)


def getindexbounds(f, eta):
    """
    Get indices of values closest below and above maximum value in an array
    :param f: array
    :type f: `~numpy.ndarray`
    :param eta: look for value in array that is closes to max_value * eta
    :type eta: float
    :return: tuple containing index of max value, index of value closest below max value,
     index of value closest above max value
    :rtype: (int, int, int)
    """
    mi = f.argmax()  # get indices of max values
    m = max(f)  # get maximum value
    b = m * eta  #
    l = find_nearest(f[:mi], b)  # find closest value below max value
    u = find_nearest(f[mi:], b) + mi  # find closest value above max value
    return mi, l, u


def gen_Pool(ncores=0):
    """
    Generate mulitprocessing pool object utilizing ncores amount of cores
    :param ncores: number of CPU cores for multiprocessing.Pool, if ncores == 0 use all available
    :type ncores: int
    :return: multiprocessing.Pool object
    :rtype: `~multiprocessing.Pool`
    """
    import multiprocessing

    ncores_max = multiprocessing.cpu_count()

    if ncores == 0 or ncores > ncores_max:
        ncores = ncores_max
    if ncores > ncores_max:
        print('Reduced number of requested CPU slots to available number: {}'.format(ncores))

    print('gen_Pool: Generated multiprocessing Pool with {} cores\n'.format(ncores))

    pool = multiprocessing.Pool(ncores)
    return pool


def excludeQualityClasses(picks, qClasses, timeerrorsP, timeerrorsS):
    """
    takes PyLoT picks dictionary and returns a new dictionary with certain classes excluded.
    :param picks: PyLoT picks dictionary
    :type picks: dict
    :param qClasses: list (or int) of quality classes (0-4) to exclude
    :type qClasses: [int]
    :param timeerrorsP: width of quality classes for P onsets in seconds
    :type timeerrorsP: (float, float, float, float)
    :param timeerrorsS: width of quality classes for S onsets in seconds
    :type timeerrorsS: (float, float, float, float])
    :return: dictionary containing only picks above the excluded quality class(es)
    :rtype: dict
    """
    from pylot.core.pick.utils import get_quality_class

    if type(qClasses) in [int, float]:
        qClasses = [qClasses]

    picksdict_new = {}

    phaseError = {'P': timeerrorsP,
                  'S': timeerrorsS}

    for station, phases in picks.items():
        for phase, pick in phases.items():
            if not type(pick) in [AttribDict, dict]:
                continue
            pickerror = phaseError[identifyPhaseID(phase)]
            quality = get_quality_class(pick['spe'], pickerror)
            if not quality in qClasses:
                if not station in picksdict_new:
                    picksdict_new[station] = {}
                picksdict_new[station][phase] = pick

    return picksdict_new


def clims(lim1, lim2):
    """
    takes two pairs of limits and returns one pair of common limts
    :param lim1: limit 1
    :type lim1: int
    :param lim2: limit 2
    :type lim2: int
    :return: new upper and lower limit common to both given limits
    :rtype: [int, int]

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
    :param window: time window whitin which data is demeaned
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


def find_in_list(list, str):
    """
    takes a list of strings and a string and returns the first list item
    matching the string pattern
    :param list: list to search in
    :type list: list
    :param str: pattern to search for
    :type str: str
    :return: first list item containing pattern
    :rtype: str

    .. example::

    >>> l = ['/dir/e1234.123.12', '/dir/e2345.123.12', 'abc123', 'def456']
    >>> find_in_list(l, 'dir')
    '/dir/e1234.123.12'
    >>> find_in_list(l, 'e1234')
    '/dir/e1234.123.12'
    >>> find_in_list(l, 'e2')
    '/dir/e2345.123.12'
    >>> find_in_list(l, 'ABC')
    'abc123'
    >>> find_in_list(l, 'f456')
    'def456'
    >>> find_in_list(l, 'gurke')

    """
    rlist = [s for s in list if str.lower() in s.lower()]
    if rlist:
        return rlist[0]
    return None


def find_nearest(array, value):
    """
    function find_nearest takes an array and a value and returns the
    index of the nearest value found in the array
    :param array: array containing values
    :type array: `~numpy.ndarray`
    :param value: number searched for
    :type value: float
    :return: index of the array item being nearest to the value
    :rtype: int

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
    """
    return (np.abs(array - value)).argmin()


def fnConstructor(s):
    """
    takes a string and returns a valid filename (especially on windows machines)
    :param s: desired filename
    :type s: str
    :return: valid filename
    :rtype: str
    """
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


def get_None(value):
    """
    Convert "None" to None
    :param value:
    :type value: str, bool
    :return:
    :rtype: bool
    """
    if value == 'None':
        return None
    else:
        return value


def get_Bool(value):
    """
    Convert string representations of bools to their true boolean value
    :param value:
    :type value: str, bool
    :return: true boolean value
    :rtype: bool
    """
    if value in ['True', 'true']:
        return True
    elif value in ['False', 'false']:
        return False
    else:
        return value


def four_digits(year):
    """
    takes a two digit year integer and returns the correct four digit equivalent
    from the last 100 years
    :param year: two digit year
    :type year: int
    :return: four digit year correspondent
    :rtype: int

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


def common_range(stream):
    """
    takes a stream object and returns the earliest end and the latest start time of all contained trace objects
    :param stream:  seismological data stream
    :type stream: `~obspy.core.stream.Stream`
    :return: maximum start time and minimum end time
    :rtype: (`~maximum start time and minimum end time`, maximum start time and minimum end time)
    """
    max_start = None
    min_end = None
    for trace in stream:
        if max_start is None or trace.stats.starttime > max_start:
            max_start = trace.stats.starttime
        if min_end is None or trace.stats.endtime < min_end:
            min_end = trace.stats.endtime
    return max_start, min_end


def full_range(stream):
    """
    takes a stream object and returns the latest end and the earliest start
    time of all contained trace objects
    :param stream: seismological data stream
    :type stream: `~obspy.core.stream.Stream`
    :return: minimum start time and maximum end time
    :rtype: (`~maximum start time and minimum end time`, maximum start time and minimum end time)
    """
    if not stream:
        print('full_range: Empty Stream!')
        return None, None

    min_start = min([trace.stats.starttime for trace in stream])
    max_end = max([trace.stats.endtime for trace in stream])

    return min_start, max_end


def transformFilteroptions2String(filtopts):
    st = ''
    if not filtopts:
        return st
    if 'type' in filtopts.keys():
        st += '{}'.format(filtopts['type'])
        if 'freq' in filtopts.keys():
            st += ' | freq: {}'.format(filtopts['freq'])
        elif 'freqmin' in filtopts.keys() and 'freqmax' in filtopts.keys():
            st += ' | freqmin: {} | freqmax: {}'.format(filtopts['freqmin'], filtopts['freqmax'])
    for key, value in filtopts.items():
        if key in ['type', 'freq', 'freqmin', 'freqmax']:
            continue
        st += ' | {}: {}'.format(key, value)
    return st


def transformFilterString4Export(st):
    st = st.replace('|', '//')
    st = st.replace(':', '/')
    st = st.replace(' ', '')
    return st


def backtransformFilterString(st):
    st = st.split('smi:local/')
    st = st[1] if len(st) > 1 else st[0]
    st = st.replace('//', ' | ')
    st = st.replace('/', ': ')
    return st


def getHash(time):
    """
    takes a time object and returns the corresponding SHA1 hash of the formatted date string
    :param time: time object for which a hash should be calculated
    :type time: `~obspy.core.utcdatetime.UTCDateTime`
    :return: SHA1 hash
    :rtype: str
    """
    hg = hashlib.sha1()
    hg.update(time.strftime('%Y-%m-%d %H:%M:%S.%f'))
    return hg.hexdigest()


def getLogin():
    """
    returns the actual user's login ID
    :return: login ID
    :rtype: str
    """
    import getpass
    return getpass.getuser()


def getOwner(fn):
    """
    takes a filename and return the login ID of the actual owner of the file
    :param fn: filename of the file tested
    :type fn: str
    :return: login ID of the file's owner
    :rtype: str
    """
    system_name = platform.system()
    if system_name in ["Linux", "Darwin"]:
        import pwd
        return pwd.getpwuid(os.stat(fn).st_uid).pw_name
    elif system_name == "Windows":
        import win32security
        f = win32security.GetFileSecurity(fn, win32security.OWNER_SECURITY_INFORMATION)
        (username, domain, sid_name_use) = win32security.LookupAccountSid(None, f.GetSecurityDescriptorOwner())
        return username


def getPatternLine(fn, pattern):
    """
    takes a file name and a pattern string to search for in the file and
    returns the first line which contains the pattern string otherwise 'None'
    :param fn: file name
    :type fn: str
    :param pattern: pattern string to search for
    :type pattern: str
    :return: the complete line containing the pattern string or None
    :rtype: int, None


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


def is_executable(fn):
    """
    takes a filename and returns True if the file is executable on the system
    and False otherwise
    :param fn: path to the file to be tested
    :return: True or False
    :rtype: bool
    """
    return os.path.isfile(fn) and os.access(fn, os.X_OK)


def isSorted(iterable):
    """
    takes an iterable and returns True if the items are in order otherwise False
    :param iterable: an iterable object
    :type iterable:
    :return: Boolean
    :rtype: bool

    ..example::
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
    """
    assert isIterable(iterable), 'object is not iterable; object: {' \
                                 '}'.format(iterable)
    if type(iterable) is str:
        iterable = [s for s in iterable]
    return sorted(iterable) == iterable


def isIterable(obj):
    """
    takes a python object and returns True is the object is iterable and
    False otherwise
    :param obj: a python object
    :type obj: object
    :return: True of False
    :rtype: bool
    """
    try:
        iterator = iter(obj)
    except TypeError as te:
        return False
    return True


def key_for_set_value(d):
    """
    takes a dictionary and returns the first key for which's value the
    boolean is True
    :param d: dictionary containing values
    :type d: dict
    :return: key to the first non-False value found; None if no value's
    boolean equals True
    :rtype:
    """
    r = None
    for k, v in d.items():
        if v:
            return k
    return r


def prepTimeAxis(stime, trace, verbosity=0):
    """
    takes a starttime and a trace object and returns a valid time axis for
    plotting
    :param stime: start time of the actual seismogram as UTCDateTime
    :type stime: `~obspy.core.utcdatetime.UTCDateTime`
    :param trace: seismic trace object
    :type trace: `~obspy.core.trace.Trace`
    :param verbosity: if != 0, debug output will be written to console
    :type verbosity: int
    :return: valid numpy array with time stamps for plotting
    :rtype: `~numpy.ndarray`
    """
    nsamp = trace.stats.npts
    srate = trace.stats.sampling_rate
    tincr = trace.stats.delta
    etime = stime + nsamp / srate
    time_ax = np.linspace(stime, etime, nsamp)
    if len(time_ax) < nsamp:
        if verbosity:
            print('elongate time axes by one datum')
        time_ax = np.arange(stime, etime + tincr, tincr)
    elif len(time_ax) > nsamp:
        if verbosity:
            print('shorten time axes by one datum')
        time_ax = np.arange(stime, etime - tincr, tincr)
    if len(time_ax) != nsamp:
        print('Station {0}, {1} samples of data \n '
              '{2} length of time vector \n'
              'delta: {3}'.format(trace.stats.station,
                                  nsamp, len(time_ax), tincr))
        time_ax = None
    return time_ax


def find_horizontals(data):
    """
    takes `obspy.core.stream.Stream` object and returns a list containing the component labels of the horizontal components available
    :param data: waveform data
    :type data: `obspy.core.stream.Stream`
    :return: components list
    :rtype: list

    ..example::

    >>> st = read()
    >>> find_horizontals(st)
    [u'N', u'E']
    """
    rval = []
    for tr in data:
        if tr.stats.channel[-1].upper() in ['Z', '3']:
            continue
        else:
            rval.append(tr.stats.channel[-1].upper())
    return rval


def pick_color(picktype, phase, quality=0):
    """
    Create pick color by modifying the base color by the quality.

    Returns rgba values in a range of [0, 255]. picktype, phase decide the base color,
    quality decides the applied modifier
    :param picktype: 'manual' or 'automatic'
    :type picktype: str
    :param phase: 'P' or 'S'
    :type phase: str
    :param quality: quality of pick. Decides the new intensity of the modifier color
    :type quality: int
    :return: tuple containing modified rgba color values
    :rtype: (int, int, int, int)
    """
    min_quality = 3
    bpc = base_phase_colors(picktype, phase)  # returns dict like {'modifier': 'g', 'rgba': (0, 0, 255, 255)}
    rgba = bpc['rgba']
    modifier = bpc['modifier']
    intensity = 255. * quality / min_quality
    rgba = modify_rgba(rgba, modifier, intensity)
    return rgba


def pick_color_plt(picktype, phase, quality=0):
    """
    Create pick color by modifying the base color by the quality.

    Returns rgba values in a range of [0, 1]. picktype, phase decide the base color,
    quality decides the applied modifier
    :param picktype: 'manual' or 'automatic'
    :type picktype: str
    :param phase: 'P' or 'S'
    :type phase: str
    :param quality: quality of pick. Decides the new intensity of the modifier color
    :type quality: int
    :return: tuple containing rgba values matplotlib style, ranging from [0, 1]
    :rtype: (float, float, float, float)
    """
    rgba = list(pick_color(picktype, phase, quality))
    for index, val in enumerate(rgba):
        rgba[index] /= 255.
    return rgba


def pick_linestyle_plt(picktype, key):
    """
    Get matplotlib line style for plotting by pick type and pick parameter (earliest/latest possible pick,
    symmetric picking error or most probable pick).
    :param picktype: 'manual' or 'automatic'
    :type picktype: str
    :param key: which pick parameter should be plotted, 'mpp', 'epp', 'lpp' or 'spe'
    :type key: str
    :return: tuple containing matplotlib line style string and line thicknes
    :rtype: (str, float)
    """
    linestyles_manu = {'mpp': ('solid', 2.),
                       'epp': ('dashed', 1.),
                       'lpp': ('dashed', 1.),
                       'spe': ('dashed', 1.)}
    linestyles_auto = {'mpp': ('dotted', 2.),
                       'epp': ('dashdot', 1.),
                       'lpp': ('dashdot', 1.),
                       'spe': ('dashdot', 1.)}
    linestyles = {'manual': linestyles_manu,
                  'auto': linestyles_auto}
    return linestyles[picktype][key]


def modify_rgba(rgba, modifier, intensity):
    """
    Modify rgba color by adding the given intensity to the modifier color
    :param rgba: tuple containing rgba values
    :type rgba: (int, int, int, int)
    :param modifier: which color should be modified, eg. 'r', 'g', 'b'
    :type modifier: str
    :param intensity: intensity to be added to selected color
    :type intensity: float
    :return: tuple containing rgba values
    :rtype: (int, int, int, int)
    """
    rgba = list(rgba)
    index = {'r': 0,
             'g': 1,
             'b': 2}
    val = rgba[index[modifier]] + intensity
    if val > 255.:
        val = 255.
    elif val < 0.:
        val = 0
    rgba[index[modifier]] = val
    return tuple(rgba)


def base_phase_colors(picktype, phase):
    """
    Get base color for a phase from style settings
    :param picktype: 'manual' or 'automatic' picks
    :type picktype: str
    :param phase: Phase to select color for, 'P' or 'S'
    :type phase: str
    :return: dictionary {'modifier': 'g', 'rgba': (0, 0, 255, 255)}
    :rtype: dict
    """
    phasecolors = style_settings.phasecolors
    return phasecolors[picktype][phase]


def transform_colors_mpl_str(colors, no_alpha=False):
    """
    Transforms rgba color values to a matplotlib string of color values with a range of [0, 1]
    :param colors: tuple of rgba color values ranging from [0, 255]
    :type colors: (float, float, float, float)
    :param no_alpha: Wether to return a alpha value in the matplotlib color string
    :type no_alpha: bool
    :return: String containing r, g, b values and alpha value if no_alpha is False (default)
    :rtype: str
    """
    colors = list(colors)
    colors_mpl = tuple([color / 255. for color in colors])
    if no_alpha:
        colors_mpl = '({}, {}, {})'.format(*colors_mpl)
    else:
        colors_mpl = '({}, {}, {}, {})'.format(*colors_mpl)
    return colors_mpl


def transform_colors_mpl(colors):
    """
    Transform rgba colors from [0, 255] to [0, 1]
    :param colors: tuple of rgba color values ranging from [0, 255]
    :type colors: (float, float, float, float)
    :return: tuple of rgba color values ranging from [0, 1]
    :rtype: (float, float, float, float)
    """
    colors = list(colors)
    colors_mpl = tuple([color / 255. for color in colors])
    return colors_mpl


def remove_underscores(data):
    """
    takes a `obspy.core.stream.Stream` object and removes all underscores
    from station names
    :param data: stream of seismic data
    :type data: `~obspy.core.stream.Stream`
    :return: data stream
    :rtype: `~obspy.core.stream.Stream`
    """
    # for tr in data:
    #    # remove underscores
    #    tr.stats.station = tr.stats.station.strip('_')
    return data


def trim_station_components(data, trim_start=True, trim_end=True):
    """
    cut a stream so only the part common to all three traces is kept to avoid dealing with offsets
    :param data: stream of seismic data
    :type data: `~obspy.core.stream.Stream`
    :param trim_start: trim start of stream
    :type trim_start: bool
    :param trim_end: trim end of stream
    :type trim_end: bool
    :return: data stream
    :rtype: `~obspy.core.stream.Stream`
    """
    starttime = {False: None}
    endtime = {False: None}

    stations = get_stations(data)

    print('trim_station_components: Will trim stream for trim_start: {} and for '
          'trim_end: {}.'.format(trim_start, trim_end))
    for station in stations:
        wf_station = data.select(station=station)
        starttime[True] = max([trace.stats.starttime for trace in wf_station])
        endtime[True] = min([trace.stats.endtime for trace in wf_station])
        wf_station.trim(starttime=starttime[trim_start], endtime=endtime[trim_end])

    return data


def merge_stream(stream):
    gaps = stream.get_gaps()
    if gaps:
        # list of merged stations (seed_ids)
        merged = ['{}.{}.{}.{}'.format(*gap[:4]) for gap in gaps]
        stream.merge(method=1)
        print('Merged the following stations because of gaps:')
        for merged_station in merged:
            print(merged_station)

    return stream, gaps


def check4gapsAndRemove(data):
    """
    check for gaps in Stream and remove them
    :param data: stream of seismic data
    :type data: `~obspy.core.stream.Stream`
    :return: data stream
    :rtype: `~obspy.core.stream.Stream`
    """
    stations = get_stations(data)

    for station in stations:
        wf_station = data.select(station=station)
        if wf_station.get_gaps():
            for trace in wf_station:
                data.remove(trace)
            print('check4gaps: Found gaps and removed station {} from waveform data.'.format(station))

    return data


def check4gapsAndMerge(data):
    """
    check for gaps in Stream and merge if gaps are found
    :param data: stream of seismic data
    :type data: `~obspy.core.stream.Stream`
    :return: data stream
    :rtype: `~obspy.core.stream.Stream`
    """
    gaps = data.get_gaps()
    if gaps:
        merged = ['{}.{}.{}.{}'.format(*gap[:4]) for gap in gaps]
        data.merge(method=1)
        print('Merged the following stations because of gaps:')
        for merged_station in merged:
            print(merged_station)

    return data


def check4doubled(data):
    """
    check for doubled stations for same channel in Stream and take only the first one
    :param data: stream of seismic data
    :type data: `~obspy.core.stream.Stream`
    :return: data stream
    :rtype: `~obspy.core.stream.Stream`
    """
    stations = get_stations(data)

    for station in stations:
        wf_station = data.select(station=station)
        # create list of all possible channels
        channels = []
        for trace in wf_station:
            channel = trace.stats.channel
            if not channel in channels:
                channels.append(channel)
            else:
                print('check4doubled: removed the following trace for station {}, as there is'
                      ' already a trace with the same channel given:\n{}'.format(
                    station, trace
                ))
                data.remove(trace)
    return data


def get_stations(data):
    """
    Get list of all station names in data stream
    :param data: stream containing seismic traces
    :type data: `~obspy.core.stream.Stream`
    :return: list of all station names in data, no duplicates
    :rtype: list of str
    """
    stations = []
    for tr in data:
        station = tr.stats.station
        if not station in stations:
            stations.append(station)

    return stations


def check4rotated(data, metadata=None, verbosity=1):
    """
    Check all traces in data. If a trace is not in ZNE rotation (last symbol of channel code is numeric) and the trace
    is in the metadata with azimuth and dip, rotate it to classical ZNE orientation.
    Rotating the traces requires them to be of the same length, so, all traces will be trimmed to a common length as a
    side effect.
    :param data: stream object containing seismic traces
    :type data: `~obspy.core.stream.Stream`
    :param metadata: tuple containing metadata type string and metadata parser object
    :type metadata: (str, `~obspy.io.xseed.parser.Parser`)
    :param verbosity: if 0 print no information at runtime
    :type verbosity: int
    :return: stream object with traditionally oriented traces (ZNE) for stations that had misaligned traces (123) before
    :rtype: `~obspy.core.stream.Stream`
    """

    def rotate_components(wfstream, metadata=None):
        """
        Rotate components if orientation code is numeric (= non traditional orientation).

        Azimut and dip are fetched from metadata. To be rotated, traces of a station have to be cut to the same length.
        Returns unrotated traces of no metadata is provided
        :param wfstream: stream containing seismic traces of a station
        :type wfstream: `~obspy.core.stream.Stream`
        :param metadata: tuple containing metadata type string and metadata parser object
        :type metadata: (str, `~obspy.io.xseed.parser.Parser`)
        :return: stream object with traditionally oriented traces (ZNE)
        :rtype: `~obspy.core.stream.Stream`
        """

        # check if any traces in this station need to be rotated
        trace_ids = [trace.id for trace in wfstream]
        orientations = [trace_id[-1] for trace_id in trace_ids]
        rotation_required = [orientation.isnumeric() for orientation in orientations]
        if any(rotation_required):
            t_start = full_range(wfstream)
            try:
                azimuts = [metadata.get_coordinates(tr_id, t_start)['azimuth'] for tr_id in trace_ids]
                dips = [metadata.get_coordinates(tr_id, t_start)['dip'] for tr_id in trace_ids]
            except (KeyError, TypeError) as e:
                print('Failed to rotate trace {}, no azimuth or dip available in metadata'.format(trace_id))
                return wfstream
            if len(wfstream) < 3:
                print('Failed to rotate Stream {}, not enough components available.'.format(wfstream))
                return wfstream
            # to rotate all traces must have same length, so trim them
            wfstream = trim_station_components(wfstream, trim_start=True, trim_end=True)
            try:
                z, n, e = rotate2zne(wfstream[0], azimuts[0], dips[0],
                                     wfstream[1], azimuts[1], dips[1],
                                     wfstream[2], azimuts[2], dips[2])
                print('check4rotated: rotated trace {} to ZNE'.format(trace_id))
                # replace old data with rotated data, change the channel code to ZNE
                z_index = dips.index(min(
                    dips))  # get z-trace index, z has minimum dip of -90 (dip is measured from 0 to -90, with -90 being vertical)
                wfstream[z_index].data = z
                wfstream[z_index].stats.channel = wfstream[z_index].stats.channel[0:-1] + 'Z'
                del trace_ids[z_index]
                for trace_id in trace_ids:
                    coordinates = metadata.get_coordinates(trace_id, t_start)
                    dip, az = coordinates['dip'], coordinates['azimuth']
                    trace = wfstream.select(id=trace_id)[0]
                    if az > 315 or az <= 45 or az > 135 and az <= 225:
                        trace.data = n
                        trace.stats.channel = trace.stats.channel[0:-1] + 'N'
                    elif az > 45 and az <= 135 or az > 225 and az <= 315:
                        trace.data = e
                        trace.stats.channel = trace.stats.channel[0:-1] + 'E'
            except (ValueError) as e:
                print(e)
                return wfstream

        return wfstream

    if metadata is None:
        if verbosity:
            msg = 'Warning: could not rotate traces since no metadata was given\nset Inventory file!'
            print(msg)
        return data
    stations = get_stations(data)
    for station in stations:  # loop through all stations and rotate data if neccessary
        wf_station = data.select(station=station)
        rotate_components(wf_station, metadata)
    return data


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

    subprocess.check_output('{} | tee /dev/stderr'.format(cmd), shell=True)


def loopIdentifyPhase(phase):
    """
    Loop through phase string and try to recognize its type (P or S wave).
    Global variable ALTSUFFIX gives alternative suffix for phases if they do not end with P, p or S, s.
    If ALTSUFFIX is not given, the function will cut the last letter of the phase string until string ends
    with P or S.
    :param phase: phase name
    :type phase: str
    :return: str of phase ending with identified type, None if phase could not be identified
    :rtype: str or None
    """
    from pylot.core.util.defaults import ALTSUFFIX

    if phase is None:
        raise NameError('Can not identify phase that is None')

    phase_copy = phase
    while not identifyPhase(phase_copy):
        identified = False
        for alt_suf in ALTSUFFIX:
            if phase_copy.endswith(alt_suf):
                phase_copy = phase_copy.split(alt_suf)[0]
                identified = True
        if not identified:
            phase_copy = phase_copy[:-1]
        if len(phase_copy) < 1:
            print('Warning: Could not identify phase {}!'.format(phase))
            return
    return phase_copy


def identifyPhase(phase):
    """
    Returns capital P or S if phase string is identified by last letter. Else returns False.
    :param phase: phase name
    :type phase: str
    :return: 'P', 'S' or False
    :rtype: str or bool
    """
    # common phase suffix for P and S
    common_P = ['P', 'p', 'R']
    common_S = ['S', 's']
    if phase is None:
        return False
    if phase[-1] in common_P:
        return 'P'
    if phase[-1] in common_S:
        return 'S'
    else:
        return False


def identifyPhaseID(phase):
    """
    Returns phase id (capital P or S)
    :param phase: phase name
    :type phase: str
    :return: phase type string
    :rtype: str
    """
    return identifyPhase(loopIdentifyPhase(phase))


def check_all_obspy(eventlist):
    ev_type = 'obspydmt'
    return check_event_folders(eventlist, ev_type)


def check_all_pylot(eventlist):
    ev_type = 'pylot'
    return check_event_folders(eventlist, ev_type)


def check_event_folders(eventlist, ev_type):
    checklist = []
    clean_eventlist = []
    for path in eventlist:
        folder_check = check_event_folder(path)
        if not folder_check:
            warnings.warn('Unrecognized event folder: {}'.format(path))
            continue
        checklist.append(folder_check == ev_type)
        clean_eventlist.append(path)
    if all(checklist) or len(checklist) == 0:
        return clean_eventlist
    else:
        warnings.warn('Not all selected folders of type {}'.format(ev_type))
        return []


def check_event_folder(path):
    ev_type = None
    folder = path.split('/')[-1]
    # for pylot: select only folders that start with 'e', containin two dots and have length 12
    if (folder.startswith('e')
            and len(folder.split('.')) == 3
            and len(folder) == 12):
        ev_type = 'pylot'
    elif check_obspydmt_eventfolder(folder)[0]:
        ev_type = 'obspydmt'
    return ev_type


def correct_iplot(iplot):
    """
    iplot should be in range 0...2, but it can be given as True or 'True' as well, which should be converted
    to an integer. Both will be converted to 2.
    :type iplot: Bool or int
    :return: iplot as an integer
    :rtype: int
    """
    # TODO this is a hack, there should never be the ability to pass anything else but an int
    try:
        iplot = int(iplot)
    except ValueError:
        if get_Bool(iplot):
            iplot = 2
        else:
            iplot = 0
    return iplot


def station_id_remove_channel(station_id):
    """
    Remove the channel from a SEED station id and return Network.Station.Location.
    >>> station_id_remove_channel("BW.MANZ..EHZ")
    'BW.MANZ.'
    >>> station_id_remove_channel("BW.MANZ.A.EHZ")
    'BW.MANZ.A'

    :param station_id:
    :return: station id with channel removed
    """
    # split at the last occuring dot and keep the left part
    station_id = station_id.rpartition('.')[0]
    return station_id


class SetChannelComponents(object):
    def __init__(self):
        self.setDefaultCompPosition()

    def setDefaultCompPosition(self):
        # default component order
        self.compPosition_Map = dict(Z=2, N=1, E=0)
        self.compName_Map = {'3': 'Z',
                             '1': 'N',
                             '2': 'E'}

    def _getCurrentPosition(self, component):
        for key, value in self.compName_Map.items():
            if value == component:
                return key, value
        errMsg = 'getCurrentPosition: Could not find former position of component {}.'.format(component)
        raise ValueError(errMsg)

    def _switch(self, component, component_alter, verbosity=0):
        # Without switching, multiple definitions of the same alter_comp are possible
        old_alter_comp, _ = self._getCurrentPosition(component)
        old_comp = self.compName_Map[component_alter]
        if not old_alter_comp == component_alter and not old_comp == component:
            self.compName_Map[old_alter_comp] = old_comp
            if verbosity > 0:
                print('switch: Automatically switched component {} to {}'.format(old_alter_comp, old_comp))

    def setCompPosition(self, component_alter, component, switch=True, verbosity=1):
        component_alter = str(component_alter)
        if not component_alter in self.compName_Map.keys():
            errMsg = 'setCompPosition: Unrecognized alternative component {}. Expecting one of {}.'
            raise ValueError(errMsg.format(component_alter, self.compName_Map.keys()))
        if not component in self.compPosition_Map.keys():
            errMsg = 'setCompPosition: Unrecognized target component {}. Expecting one of {}.'
            raise ValueError(errMsg.format(component, self.compPosition_Map.keys()))
        if verbosity > 0:
            print('setCompPosition: set component {} to {}'.format(component_alter, component))
        if switch:
            self._switch(component, component_alter, verbosity)
        self.compName_Map[component_alter] = component

    def getCompPosition(self, component):
        return self._getCurrentPosition(component)[0]

    def getPlotPosition(self, component):
        component = str(component)
        if component in self.compPosition_Map.keys():
            return self.compPosition_Map[component]
        elif component in self.compName_Map.keys():
            return self.compPosition_Map[self.compName_Map[component]]
        else:
            errMsg = 'getCompPosition: Unrecognized component {}. Expecting one of {} or {}.'
            raise ValueError(errMsg.format(component, self.compPosition_Map.keys(), self.compName_Map.keys()))

    @staticmethod
    def from_qsettings(settings):
        scc = SetChannelComponents()
        for value in ['Z', 'N', 'E']:
            key = settings.value(value, None)
            if not key:
                print('Could not get channel component map from QSettings. Writing default channel order to QSettings.')
                scc.setDefaultCompPosition()
                for value in ['Z', 'N', 'E']:
                    settings.setValue(value, scc.getCompPosition(value))
                return scc
            scc.setCompPosition(key, value, verbosity=0)
        return scc


if __name__ == "__main__":
    import doctest

    doctest.testmod()
