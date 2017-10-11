#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hashlib
import os
import platform
import re
import subprocess

import numpy as np
from obspy import UTCDateTime, read
from obspy.core import AttribDict
from obspy.signal.rotate import rotate2zne
from obspy.io.xseed.utils import SEEDParserException

from pylot.core.io.inputs import PylotParameter
from pylot.styles import style_settings

from scipy.interpolate import splrep, splev
from PySide import QtCore, QtGui

try:
    import pyqtgraph as pg
except Exception as e:
    print('PyLoT: Could not import pyqtgraph. {}'.format(e))
    pg = None

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)


def readDefaultFilterInformation(fname):
    pparam = PylotParameter(fname)
    return readFilterInformation(pparam)


def readFilterInformation(pylot_parameter):
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
    return splev, splrep(x, y)


def getindexbounds(f, eta):
    mi = f.argmax()
    m = max(f)
    b = m * eta
    l = find_nearest(f[:mi], b)
    u = find_nearest(f[mi:], b) + mi
    return mi, l, u


def gen_Pool(ncores=0):
    '''
    :param ncores: number of CPU cores for multiprocessing.Pool, if ncores == 0 use all available
    :return: multiprocessing.Pool object
    '''
    import multiprocessing

    if ncores == 0:
        ncores = multiprocessing.cpu_count()

    print('gen_Pool: Generated multiprocessing Pool with {} cores\n'.format(ncores))

    pool = multiprocessing.Pool(ncores)
    return pool


def excludeQualityClasses(picks, qClasses, timeerrorsP, timeerrorsS):
    '''
    takes PyLoT picks dictionary and returns a new dictionary with certain classes excluded.
    :param picks: PyLoT picks dictionary
    :param qClasses: list (or int) of quality classes (0-4) to exclude
    :param timeerrorsP: time errors for classes (0-4) for P
    :param timeerrorsS: time errors for classes (0-4) for S
    :return: new picks dictionary
    '''
    from pylot.core.pick.utils import getQualityFromUncertainty

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
            quality = getQualityFromUncertainty(pick['spe'], pickerror)
            if not quality in qClasses:
                if not station in picksdict_new:
                    picksdict_new[station] = {}
                picksdict_new[station][phase] = pick

    return picksdict_new


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


def find_in_list(list, str):
    """
    takes a list of strings and a string and returns the first list item
    matching the string pattern
    :param list: list to search in
    :param str: pattern to search for
    :return: first list item containing pattern

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


def real_None(value):
    if value == 'None':
        return None
    else:
        return value


def real_Bool(value):
    if value == 'True':
        return True
    elif value == 'False':
        return False
    else:
        return value


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


def common_range(stream):
    '''
    takes a stream object and returns the earliest end and the latest start
    time of all contained trace objects
    :param stream: seismological data stream
    :type stream: `~obspy.core.stream.Stream`
    :return: maximum start time and minimum end time
    '''
    max_start = None
    min_end = None
    for trace in stream:
        if max_start is None or trace.stats.starttime > max_start:
            max_start = trace.stats.starttime
        if min_end is None or trace.stats.endtime < min_end:
            min_end = trace.stats.endtime
    return max_start, min_end


def full_range(stream):
    '''
    takes a stream object and returns the latest end and the earliest start
    time of all contained trace objects
    :param stream: seismological data stream
    :type stream: `~obspy.core.stream.Stream`
    :return: minimum start time and maximum end time
    '''
    min_start = min([trace.stats.starttime for trace in stream])
    max_end = max([trace.stats.endtime for trace in stream])

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
    return os.getlogin()


def getOwner(fn):
    '''
    takes a filename and return the login ID of the actual owner of the file
    :param fn: filename of the file tested
    :type fn: str
    :return: login ID of the file's owner
    '''
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
    """
    return os.path.isfile(fn) and os.access(fn, os.X_OK)


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


def key_for_set_value(d):
    """
    takes a dictionary and returns the first key for which's value the
    boolean is True
    :param d: dictionary containing values
    :type d: dict
    :return: key to the first non-False value found; None if no value's
    boolean equals True
    """
    r = None
    for k, v in d.items():
        if v:
            return k
    return r


def prepTimeAxis(stime, trace, verbosity=0):
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


def make_pen(picktype, phase, key, quality):
    if pg:
        rgba = pick_color(picktype, phase, quality)
        linestyle, width = pick_linestyle_pg(picktype, key)
        pen = pg.mkPen(rgba, width=width, style=linestyle)
        return pen


def pick_color(picktype, phase, quality=0):
    min_quality = 3
    bpc = base_phase_colors(picktype, phase)
    rgba = bpc['rgba']
    modifier = bpc['modifier']
    intensity = 255.*quality/min_quality
    rgba = modify_rgba(rgba, modifier, intensity)
    return rgba


def pick_color_plt(picktype, phase, quality=0):
    rgba = list(pick_color(picktype, phase, quality))
    for index, val in enumerate(rgba):
        rgba[index] /= 255.
    return rgba


def pick_linestyle_plt(picktype, key):
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


def pick_linestyle_pg(picktype, key):
    linestyles_manu = {'mpp': (QtCore.Qt.SolidLine, 2.),
                       'epp': (QtCore.Qt.DashLine, 1.),
                       'lpp': (QtCore.Qt.DashLine, 1.),
                       'spe': (QtCore.Qt.DashLine, 1.)}
    linestyles_auto = {'mpp': (QtCore.Qt.DotLine, 2.),
                       'epp': (QtCore.Qt.DashDotLine, 1.),
                       'lpp': (QtCore.Qt.DashDotLine, 1.),
                       'spe': (QtCore.Qt.DashDotLine, 1.)}
    linestyles = {'manual': linestyles_manu,
                  'auto': linestyles_auto}
    return linestyles[picktype][key]


def modify_rgba(rgba, modifier, intensity):
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
    phasecolors = style_settings.phasecolors
    return phasecolors[picktype][phase]

def transform_colors_mpl_str(colors, no_alpha=False):
    colors = list(colors)
    colors_mpl = tuple([color / 255. for color in colors])
    if no_alpha:
        colors_mpl = '({}, {}, {})'.format(*colors_mpl)
    else:
        colors_mpl = '({}, {}, {}, {})'.format(*colors_mpl)
    return colors_mpl

def transform_colors_mpl(colors):
    colors = list(colors)
    colors_mpl = tuple([color / 255. for color in colors])
    return colors_mpl

def remove_underscores(data):
    """
    takes a `obspy.core.stream.Stream` object and removes all underscores
    from stationnames
    :param data: stream of seismic data
    :type data: `obspy.core.stream.Stream`
    :return: data stream
    """
    for tr in data:
        # remove underscores
        tr.stats.station = tr.stats.station.strip('_')
    return data


def trim_station_components(data, trim_start=True, trim_end=True):
    '''
    cut a stream so only the part common to all three traces is kept to avoid dealing with offsets
    :param data: stream of seismic data
    :type data: `obspy.core.stream.Stream`
    :param trim_start: trim start of stream
    :type trim_start: bool
    :param trim_end: trim end of stream
    :type trim_end: bool
    :return: data stream
    '''
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


def check4gaps(data):
    '''
    check for gaps in Stream and remove them
    :param data: stream of seismic data
    :return: data stream
    '''
    stations = get_stations(data)

    for station in stations:
        wf_station = data.select(station=station)
        if wf_station.get_gaps():
            for trace in wf_station:
                data.remove(trace)
            print('check4gaps: Found gaps and removed station {} from waveform data.'.format(station))

    return data


def check4doubled(data):
    '''
    check for doubled stations for same channel in Stream and take only the first one
    :param data: stream of seismic data
    :return: data stream
    '''
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
    stations = []
    for tr in data:
        station = tr.stats.station
        if not station in stations:
            stations.append(station)

    return stations


def check4rotated(data, metadata=None, verbosity=1):

    def rotate_components(wfstream, metadata=None):
        """rotates components if orientation code is numeric.
        azimut and dip are fetched from metadata"""
        try:
            # indexing fails if metadata is None
            metadata[0]
        except:
            if verbosity:
                msg = 'Warning: could not rotate traces since no metadata was given\nset Inventory file!'
                print(msg)
            return wfstream
        if metadata[0] is None:
            # sometimes metadata is (None, (None,))
            if verbosity:
                msg = 'Warning: could not rotate traces since no metadata was given\nCheck inventory directory!'
                print(msg)
            return wfstream
        else:
            parser = metadata[1]

        def get_dip_azimut(parser, trace_id):
            """gets azimut and dip for a trace out of the metadata parser"""
            dip = None
            azimut = None
            try:
                blockettes = parser._select(trace_id)
            except SEEDParserException as e:
                print(e)
                raise ValueError
            for blockette_ in blockettes:
                if blockette_.id != 52:
                    continue
                dip = blockette_.dip
                azimut = blockette_.azimuth
                break
            if dip is None or azimut is None:
                error_msg = 'Dip and azimuth not available for trace_id {}'.format(trace_id)
                raise ValueError(error_msg)
            return dip, azimut

        trace_ids = [trace.id for trace in wfstream]
        for trace_id in trace_ids:
            orientation = trace_id[-1]
            if orientation.isnumeric():
                # misaligned channels have a number as orientation
                azimuts = []
                dips = []
                for trace_id in trace_ids:
                    try:
                        dip, azimut = get_dip_azimut(parser, trace_id)
                    except ValueError as e:
                        print(e)
                        print('Failed to rotate station {}, no azimuth or dip available in metadata'.format(trace_id))
                        return wfstream
                    azimuts.append(azimut)
                    dips.append(dip)
                # to rotate all traces must have same length
                wfstream = trim_station_components(wfstream, trim_start=True, trim_end=True)
                z, n, e = rotate2zne(wfstream[0], azimuts[0], dips[0],
                                          wfstream[1], azimuts[1], dips[1],
                                          wfstream[2], azimuts[2], dips[2])
                print('check4rotated: rotated station {} to ZNE'.format(trace_id))
                z_index = dips.index(min(dips)) # get z-trace index (dip is measured from 0 to -90
                wfstream[z_index].data = z
                wfstream[z_index].stats.channel = wfstream[z_index].stats.channel[0:-1] + 'Z'
                del trace_ids[z_index]
                for trace_id in trace_ids:
                    dip, az = get_dip_azimut(parser, trace_id)
                    trace = wfstream.select(id=trace_id)[0]
                    if az > 315 and az <= 45 or az > 135 and az <= 225:
                        trace.data = n
                        trace.stats.channel = trace.stats.channel[0:-1] + 'N'
                    elif az > 45 and az <= 135 or az > 225 and az <= 315:
                        trace.data = e
                        trace.stats.channel = trace.stats.channel[0:-1] + 'E'
                break
            else:
                continue
        return wfstream

    stations = get_stations(data)

    for station in stations:
        wf_station = data.select(station=station)
        wf_station = rotate_components(wf_station, metadata)
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


def which(program, infile=None):
    """
    takes a program name and returns the full path to the executable or None
    modified after: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    :param program: name of the desired external program
    :return: full path of the executable file
    """
    try:
        from PySide.QtCore import QSettings
        settings = QSettings()
        for key in settings.allKeys():
            if 'binPath' in key:
                os.environ['PATH'] += ':{0}'.format(settings.value(key))
        if infile is None:
            # use default parameter-file name
            bpath = os.path.join(os.path.expanduser('~'), '.pylot', 'pylot.in')
        else:
            bpath = os.path.join(os.path.expanduser('~'), '.pylot', infile)

        if os.path.exists(bpath):
            nllocpath = ":" + PylotParameter(bpath).get('nllocbin')
            os.environ['PATH'] += nllocpath
    except ImportError as e:
        print(e.message)

    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate

    return None


def loopIdentifyPhase(phase):
    '''
    Loop through phase string and try to recognize its type (P or S wave).
    Global variable ALTSUFFIX gives alternative suffix for phases if they do not end with P, p or S, s.
    If ALTSUFFIX is not given, the function will cut the last letter of the phase string until string ends
    with P or S.
    :param phase: phase name (str)
    :return:
    '''
    from pylot.core.util.defaults import ALTSUFFIX

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
    '''
    Returns capital P or S if phase string is identified by last letter. Else returns False.
    :param phase: phase name (str)
    :return: 'P', 'S' or False
    '''
    # common phase suffix for P and S
    common_P = ['P', 'p']
    common_S = ['S', 's']
    if phase[-1] in common_P:
        return 'P'
    if phase[-1] in common_S:
        return 'S'
    else:
        return False


def identifyPhaseID(phase):
    return identifyPhase(loopIdentifyPhase(phase))


def has_spe(pick):
    if not 'spe' in pick.keys():
        return None
    else:
        return pick['spe']


if __name__ == "__main__":
    import doctest

    doctest.testmod()
