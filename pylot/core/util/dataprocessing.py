#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
import sys

import numpy as np
from obspy import UTCDateTime, read_inventory, read
from obspy.io.xseed import Parser
from pylot.core.util.utils import key_for_set_value, find_in_list, \
    remove_underscores, gen_Pool


class Metadata(object):
    def __init__(self, inventory=None):
        self.inventories = []
        # saves read metadata objects (Parser/inventory) for a filename
        self.inventory_files = {}
        # saves filenames holding metadata for a seed_id
        self.seed_ids = {}
        if inventory:
            if os.path.isdir(inventory):
                self.add_inventory(inventory)
            if os.path.isfile(inventory):
                self.add_inventory_file(inventory)


    def __str__(self):
        repr = 'PyLoT Metadata object including the following inventories:\n\n'
        ntotal = len(self.inventories)
        for index, inventory in enumerate(self.inventories):
            if index < 2 or (ntotal - index) < 3:
                repr += '{}\n'.format(inventory)
            if ntotal > 4 and int(ntotal/2) == index:
                repr += '...\n'
        if ntotal > 4:
            repr += '\nTotal of {} inventories. Use Metadata.inventories to see all.'.format(ntotal)
        return repr


    def __repr__(self):
        return self.__str__()


    def add_inventory(self, path_to_inventory):
        '''
        add paths to list of inventories

        :param path_to_inventory:
        :return:
        '''
        assert (os.path.isdir(path_to_inventory)), '{} is no directory'.format(path_to_inventory)
        if not path_to_inventory in self.inventories:
            self.inventories.append(path_to_inventory)


    def add_inventory_file(self, path_to_inventory_file):
        '''
        add a single file to inventory files

        :param path_to_inventory_file:
        :return:

        '''
        assert (os.path.isfile(path_to_inventory_file)), '{} is no file'.format(path_to_inventory_file)
        self.add_inventory(os.path.split(path_to_inventory_file)[0])
        if not path_to_inventory_file in self.inventory_files.keys():
            self.read_single_file(path_to_inventory_file)


    def remove_all_inventories(self):
        self.__init__()


    def remove_inventory(self, path_to_inventory):
        '''
        remove a path from inventories list

        :param path_to_inventory:
        :return:
        '''
        if not path_to_inventory in self.inventories:
            print('Path {} not in inventories list.'.format(path_to_inventory))
            return
        self.inventories.remove(path_to_inventory)
        for filename in self.inventory_files.keys():
            if filename.startswith(path_to_inventory):
                del(self.inventory_files[filename])
        for seed_id in self.seed_ids.keys():
            if self.seed_ids[seed_id].startswith(path_to_inventory):
                del(self.seed_ids[seed_id])


    def get_metadata(self, seed_id):
        # get metadata for a specific seed_id, if not already read, try to read from inventories
        if not seed_id in self.seed_ids.keys():
            self._read_inventory_data(seed_id)
        # if seed id is not found read all inventories and try to find it there
        if not seed_id in self.seed_ids.keys():
            print('No data found for seed id {}. Trying to find it in all known inventories...'.format(seed_id))
            self.read_all()
            for inv_fname, metadata_dict in self.inventory_files.items():
                # use get_coordinates to check for seed_id
                try:
                    metadata_dict['data'].get_coordinates(seed_id)
                    self.seed_ids[seed_id] = inv_fname
                    print('Found metadata for station {}!'.format(seed_id))
                    return metadata_dict
                except Exception as e:
                    continue
            print('Could not find metadata for station {}'.format(seed_id))
            return None
        fname = self.seed_ids[seed_id]
        return self.inventory_files[fname]


    def read_all(self):
        '''
        read all metadata files found in all inventories
        :return:
        '''
        for inventory in self.inventories:
            for inv_fname in os.listdir(inventory):
                inv_fname = os.path.join(inventory, inv_fname)
                if not self.read_single_file(inv_fname):
                    continue


    def read_single_file(self, inv_fname):
        # try to read a single file as Parser/Inventory if it was not already read before
        if not inv_fname in self.inventory_files.keys():
            pass
        else:
            if not self.inventory_files[inv_fname]:
                pass
            else:
                return
        try:
            invtype, robj = self._read_metadata_file(inv_fname)
            if robj == None:
                return
        except Exception as e:
            print('Could not read file {}'.format(inv_fname))
            return
        self.inventory_files[inv_fname] = {'invtype': invtype,
                                           'data': robj}
        return True


    def get_coordinates(self, seed_id):
        metadata = self.get_metadata(seed_id)
        if not metadata:
            return
        return metadata['data'].get_coordinates(seed_id)


    def get_paz(self, seed_id, time):
        metadata = self.get_metadata(seed_id)
        if not metadata:
            return
        if metadata['invtype'] in ['dless', 'dseed']:
            return metadata['data'].get_paz(seed_id, time)
        elif metadata['invtype'] in ['resp', 'xml']:
            resp = metadata['data'].get_response(seed_id, time)
            return resp.get_paz(seed_id)


    def _read_inventory_data(self, seed_id=None):
        for inventory in self.inventories:
            if self._read_metadata_iterator(path_to_inventory=inventory, station_seed_id=seed_id):
                return


    def _read_metadata_iterator(self, path_to_inventory, station_seed_id):
        '''
        search for metadata for a specific station iteratively
        '''
        station, network, location, channel = station_seed_id.split('.')
        fnames = glob.glob(os.path.join(path_to_inventory, '*' + station_seed_id + '*'))
        if not fnames:
            # search for station name in filename
            fnames = glob.glob(os.path.join(path_to_inventory, '*' + station + '*'))
        if not fnames:
            # search for network name in filename
            fnames = glob.glob(os.path.join(path_to_inventory, '*' + network + '*'))
        if not fnames:
            print('Could not find filenames matching station name, network name or seed id')
            return
        for fname in fnames:
            if fname in self.inventory_files.keys():
                if self.inventory_files[fname]:
                    # file already read
                    continue
            invtype, robj = self._read_metadata_file(os.path.join(path_to_inventory, fname))
            try:
                robj.get_coordinates(station_seed_id)
                self.inventory_files[fname] = {'invtype': invtype,
                                               'data': robj}
                if station_seed_id in self.seed_ids.keys():
                    print('WARNING: Overwriting metadata for station {}'.format(station_seed_id))
                self.seed_ids[station_seed_id] = fname
                return True
            except Exception as e:
                continue
        print('Could not find metadata for station_seed_id {} in path {}'.format(station_seed_id, path_to_inventory))


    def _read_metadata_file(self, path_to_inventory_filename):
        '''
        function reading metadata files (either dataless seed, xml or resp)
        :param path_to_inventory_filename:
        :return: file type/ending, inventory object (Parser or Inventory)
        '''
        # functions used to read metadata for different file endings (or file types)
        read_functions = {'dless': self._read_dless,
                          'dseed': self._read_dless,
                          'xml': self._read_inventory_file,
                          'resp': self._read_inventory_file}
        file_ending = path_to_inventory_filename.split('.')[-1]
        if file_ending in read_functions.keys():
            robj, exc = read_functions[file_ending](path_to_inventory_filename)
            if exc is not None:
                raise exc
            return file_ending, robj
        # in case file endings did not match the above keys, try and error
        for file_type in ['dless', 'xml']:
            robj, exc = read_functions[file_type](path_to_inventory_filename)
            if exc is None:
                return file_type, robj
        return None, None


    def _read_dless(self, path_to_inventory):
        exc = None
        try:
            parser = Parser(path_to_inventory)
        except Exception as exc:
            parser = None
        return parser, exc


    def _read_inventory_file(self, path_to_inventory):
        exc = None
        try:
            inv = read_inventory(path_to_inventory)
        except Exception as exc:
            inv = None
        return inv, exc



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


def evt_head_check(root_dir, out_dir=None):
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
    # possible file extensions specified here:
    inv = dict(dless=dlfile, xml=invfile, resp=respfile, dseed=dlfile[:])
    if os.path.isfile(path_to_inventory):
        ext = os.path.splitext(path_to_inventory)[1].split('.')[1]
        inv[ext] += [path_to_inventory]
    else:
        for ext in inv.keys():
            inv[ext] += glob.glob1(path_to_inventory, '*.{0}'.format(ext))

    invtype = key_for_set_value(inv)

    if invtype is None:
        print("Neither dataless-SEED file, inventory-xml file nor "
              "RESP-file found!")
        print("!!WRONG CALCULATION OF SOURCE PARAMETERS!!")
        robj = None,
    elif invtype == 'dless':  # prevent multiple read of large dlsv
        print("Reading metadata information from dataless-SEED file ...")
        if len(inv[invtype]) == 1:
            fullpath_inv = os.path.join(path_to_inventory, inv[invtype][0])
            robj = Parser(fullpath_inv)
        else:
            robj = inv[invtype]
    else:
        print("Reading metadata information from inventory-xml file ...")
        robj = read_inventory(inv[invtype])
    return invtype, robj


# idea to optimize read_metadata
# def read_metadata_new(path_to_inventory):
#     metadata_objects = []
#     # read multiple files from directory
#     if os.path.isdir(path_to_inventory):
#         fnames = os.listdir(path_to_inventory)
#     # read single file
#     elif os.path.isfile(path_to_inventory):
#         fnames = [path_to_inventory]
#     else:
#         print("Neither dataless-SEED file, inventory-xml file nor "
#               "RESP-file found!")
#         print("!!WRONG CALCULATION OF SOURCE PARAMETERS!!")
#         fnames = []
#
#     for fname in fnames:
#         path_to_inventory_filename = os.path.join(path_to_inventory, fname)
#         try:
#             ftype, robj = read_metadata_file(path_to_inventory_filename)
#             metadata_objects.append((ftype, robj))
#         except Exception as e:
#             print('Could not read metadata file {} '
#                   'because of the following Exception: {}'.format(path_to_inventory_filename, e))
#     return metadata_objects



def restitute_trace(input_tuple):
    def no_metadata(tr, seed_id):
        print('no metadata file found '
              'for trace {0}'.format(seed_id))
        return tr, True

    tr, metadata, unit, force = input_tuple

    remove_trace = False

    seed_id = tr.get_id()

    mdata = metadata.get_metadata(seed_id)
    if not mdata:
        return no_metadata(tr, seed_id)

    invtype = mdata['invtype']
    inobj = mdata['data']

    # check, whether this trace has already been corrected
    if 'processing' in tr.stats.keys() \
            and np.any(['remove' in p for p in tr.stats.processing]) \
            and not force:
        print("Trace {0} has already been corrected!".format(seed_id))
        return tr, False
    stime = tr.stats.starttime
    prefilt = get_prefilt(tr)
    if invtype == 'resp':
        fresp = find_in_list(inobj, seed_id)
        if not fresp:
            return no_metadata(tr, seed_id)
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
    elif invtype == None:
        return no_metadata(tr, seed_id)
    else:
        remove_trace = True
    # apply restitution to data
    print("Correcting instrument at station %s, channel %s" \
          % (tr.stats.station, tr.stats.channel))
    try:
        if invtype in ['resp', 'dless']:
            try:
                tr.simulate(**kwargs)
            except ValueError as e:
                vmsg = '{0}'.format(e)
                print(vmsg)

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
            # restitution done to copies of data thus deleting traces
            # that failed should not be a problem
            remove_trace = True

    return tr, remove_trace


def restitute_data(data, metadata, unit='VEL', force=False, ncores=0):
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
    input_tuples = []
    for tr in data:
        input_tuples.append((tr, metadata, unit, force))
        data.remove(tr)

    pool = gen_Pool(ncores)
    result = pool.imap_unordered(restitute_trace, input_tuples)
    pool.close()

    for tr, remove_trace in result:
        if not remove_trace:
            data.traces.append(tr)

    # check if ALL traces could be restituted, take care of large datasets
    # better try restitution for smaller subsets of data (e.g. station by
    # station)

    # if len(restflag) > 0:
    #     restflag = bool(np.all(restflag))
    # else:
    #     restflag = False
    return data


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
    fc21 = fny - (fny * thi[0] / 100.)
    fc22 = fny - (fny * thi[1] / 100.)
    return (tlow[0], tlow[1], fc21, fc22)


if __name__ == "__main__":
    import doctest

    doctest.testmod()
