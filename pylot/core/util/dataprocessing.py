#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import sys

import numpy as np

from obspy import UTCDateTime, read_inventory, read
from obspy.io.xseed import Parser
from pylot.core.util.utils import key_for_set_value, find_in_list, \
    remove_underscores


def time_from_header(header):
    """
    Function takes in the second line from a .gse file and takes out the date and time from that line.
    :param header: second line from .gse file
    :type header: string
    :return: a list of integers of form [year, month, day, hour, minute, second, microsecond]
    """
    timeline = header.split(' ')
    time = timeline[1].split('/') + timeline[2].split(':')
    time = time[:-1] + time[-1].split('.')
    return [int(t) for t in time]


def check_time(datetime):
    """
    Function takes in date and time as list and validates it's values by trying to make an UTCDateTime object from it
    :param datetime: list of integers [year, month, day, hour, minute, second, microsecond]
    :type datetime: list
    :return: returns True if Values are in supposed range, returns False otherwise

    >>> check_time([1999, 01, 01, 23, 59, 59, 999000])
    True
    >>> check_time([1999, 01, 01, 23, 59, 60, 999000])
    False
    >>> check_time([1999, 01, 01, 23, 59, 59, 1000000])
    False
    >>> check_time([1999, 01, 01, 23, 60, 59, 999000])
    False
    >>> check_time([1999, 01, 01, 23, 60, 59, 999000])
    False
    >>> check_time([1999, 01, 01, 24, 59, 59, 999000])
    False
    >>> check_time([1999, 01, 31, 23, 59, 59, 999000])
    True
    >>> check_time([1999, 02, 30, 23, 59, 59, 999000])
    False
    >>> check_time([1999, 02, 29, 23, 59, 59, 999000])
    False
    >>> check_time([2000, 02, 29, 23, 59, 59, 999000])
    True
    >>> check_time([2000, 13, 29, 23, 59, 59, 999000])
    False
    """
    try:
        UTCDateTime(*datetime)
        return True
    except ValueError:
        return False


def get_file_list(root_dir):
    """
    Function uses a directorie to get all the *.gse files from it.
    :param root_dir: a directorie leading to the .gse files
    :type root_dir: string
    :return: returns a list of filenames (without path to them)
    """
    file_list = glob.glob1(root_dir, '*.gse')
    return file_list


def checks_station_second(datetime, file):
    """
    Function uses the given list to check if the parameter 'second' is set to 60 by mistake
     and sets the time correctly if so. Can only correct time if no date change would be necessary.
    :param datetime: [year, month, day, hour, minute, second, microsecond]
    :return: returns the input with the correct value for second
    """
    if datetime[5] == 60:
        if datetime[4] == 59:
            if datetime[3] == 23:
                err_msg = 'Date should be next day. ' \
                          'File not changed: {0}'.format(file)
                raise ValueError(err_msg)
            else:
                datetime[3] += 1
                datetime[4] = 0
                datetime[5] = 0
        else:
            datetime[4] += 1
            datetime[5] = 0
    return datetime


def make_time_line(line, datetime):
    """
    Function takes in the original line from a .gse file and a list of date and
    time values to make a new line with corrected date and time.
    :param line: second line from .gse file.
    :type line: string
    :param datetime: list of integers [year, month, day, hour, minute, second, microsecond]
    :type datetime: list
    :return: returns a string to write it into a file.
    """
    ins_form = '{0:02d}:{1:02d}:{2:02d}.{3:03d}'
    insertion = ins_form.format(int(datetime[3]),
                                int(datetime[4]),
                                int(datetime[5]),
                                int(datetime[6] * 1e-3))
    newline = line[:16] + insertion + line[28:]
    return newline


def evt_head_check(root_dir, out_dir = None):
    """
    A function to make sure that an arbitrary number of .gse files have correct values in their header.
    :param root_dir: a directory leading to the .gse files.
    :type root_dir: string
    :param out_dir: a directory to store the new files somwhere els.
    :return: returns nothing
    """
    if not out_dir:
        print('WARNING files are going to be overwritten!')
        inp = str(raw_input('Continue? [y/N]'))
        if not inp == 'y':
            sys.exit()
    filelist = get_file_list(root_dir)
    nfiles = 0
    for file in filelist:
        infile = open(os.path.join(root_dir, file), 'r')
        lines = infile.readlines()
        infile.close()
        datetime = time_from_header(lines[1])
        if check_time(datetime):
            continue
        else:
            nfiles += 1
            datetime = checks_station_second(datetime, file)
            print('writing ' + file)
            # write File
            lines[1] = make_time_line(lines[1], datetime)
            if not out_dir:
                out = open(os.path.join(root_dir, file), 'w')
                out.writelines(lines)
                out.close()
            else:
                out = open(os.path.join(out_dir, file), 'w')
                out.writelines(lines)
                out.close()
    print(nfiles)


def read_metadata(path_to_inventory):
    """
    take path_to_inventory and return either the corresponding list of files
    found or the Parser object for a network dataless seed volume to prevent
    read overhead for large dataless seed volumes
    :param path_to_inventory:
    :return: tuple containing a either list of files or `obspy.io.xseed.Parser`
    object and the inventory type found
    :rtype: tuple
    """
    dlfile = list()
    invfile = list()
    respfile = list()
    inv = dict(dless=dlfile, xml=invfile, resp=respfile)
    if os.path.isfile(path_to_inventory):
        ext = os.path.splitext(path_to_inventory)[1].split('.')[1]
        inv[ext] += [path_to_inventory]
    else:
        for ext in inv.keys():
            inv[ext] += glob.glob1(path_to_inventory, '*.{0}'.format(ext))

    invtype = key_for_set_value(inv)

    if invtype is None:
        raise IOError("Neither dataless-SEED file, inventory-xml file nor "
                      "RESP-file found!")
    elif invtype == 'dless':  # prevent multiple read of large dlsv
        if len(inv[invtype]) == 1:
            robj = Parser(inv[invtype][0])
        else:
            robj = inv[invtype]
    else:
        robj = inv[invtype]
    return invtype, robj


def restitute_data(data, invtype, inobj, unit='VEL', force=False):
    """
    takes a data stream and a path_to_inventory and returns the corrected
    waveform data stream
    :param data: seismic data stream
    :param invtype: type of found metadata
    :param inobj: either list of metadata files or `obspy.io.xseed.Parser`
    object
    :param unit: unit to correct for (default: 'VEL')
    :param force: force restitution for already corrected traces (default:
    False)
    :return: corrected data stream
    """

    restflag = list()

    data = remove_underscores(data)

    # loop over traces
    for tr in data:
        seed_id = tr.get_id()
        # check, whether this trace has already been corrected
        # TODO read actual value of processing key in Trace's stats instead
        # of just looking for thr key: if processing is setit doesn't
        # necessarily mean that the trace has been corrected so far but only
        # processed in some way, e.g. normalized
        if 'processing' in tr.stats.keys() \
                and np.all(['remove' in p for p in tr.stats.processing]) \
                and not force:
            print("Trace {0} has already been corrected!".format(seed_id))
            continue
        stime = tr.stats.starttime
        prefilt = get_prefilt(tr)
        if invtype == 'resp':
            fresp = find_in_list(inobj, seed_id)
            if not fresp:
                raise IOError('no response file found '
                              'for trace {0}'.format(seed_id))
            fname = fresp
            seedresp = dict(filename=fname,
                            date=stime,
                            units=unit)
            kwargs = dict(paz_remove=None, pre_filt=prefilt, seedresp=seedresp)
        elif invtype == 'dless':
            if type(inobj) is list:
                fname = Parser(find_in_list(inobj, seed_id))
            else:
                fname = inobj
            seedresp = dict(filename=fname,
                            date=stime,
                            units=unit)
            kwargs = dict(pre_filt=prefilt, seedresp=seedresp)
        elif invtype == 'xml':
            invlist = inobj
            if len(invlist) > 1:
                finv = find_in_list(invlist, seed_id)
            else:
                finv = invlist[0]
            inventory = read_inventory(finv, format='STATIONXML')
        else:
            restflag.append(False)
            continue
        # apply restitution to data
        try:
            if invtype in ['resp', 'dless']:
                tr.simulate(**kwargs)
            else:
                tr.attach_response(inventory)
                tr.remove_response(output=unit,
                                   pre_filt=prefilt)
        except ValueError as e:
            msg0 = 'Response for {0} not found in Parser'.format(seed_id)
            msg1 = 'evalresp failed to calculate response'
            if msg0 not in e.message or msg1 not in e.message:
                raise
            else:
                restflag.append(False)
                continue
        restflag.append(True)
    # check if ALL traces could be restituted, take care of large datasets
    # better try restitution for smaller subsets of data (e.g. station by
    # station)
    if len(restflag) > 0:
        restflag = bool(np.all(restflag))
    else:
        restflag = False
    return data, restflag


def get_prefilt(trace, tlow=(0.5, 0.9), thi=(5., 2.), verbosity=0):
    """
    takes a `obspy.core.stream.Trace` object, taper parameters tlow and thi and
    returns the pre-filtering corner frequencies for the cosine taper for
    further processing
    :param trace: seismic data trace
    :type trace: `obspy.core.stream.Trace`
    :param tlow: tuple or list containing the desired lower corner
    frequenices for a cosine taper
    :type tlow: tuple or list
    :param thi: tuple or list containing the percentage values of the
    Nyquist frequency for the desired upper corner frequencies of the cosine
    taper
    :type thi: tuple or list
    :param verbosity: verbosity level
    :type verbosity: int
    :return: pre-filt cosine taper corner frequencies
    :rtype: tuple

    ..example::

    >>> st = read()
    >>> get_prefilt(st[0])
    (0.5, 0.9, 47.5, 49.0)
    """
    if verbosity:
        print("Calculating pre-filter values for %s, %s ..." % (
            trace.stats.station, trace.stats.channel))
    # get corner frequencies for pre-filtering
    fny = trace.stats.sampling_rate / 2
    fc21 = fny - (fny * thi[0]/100.)
    fc22 = fny - (fny * thi[1]/100.)
    return (tlow[0], tlow[1], fc21, fc22)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
