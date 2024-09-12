#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import json

from obspy import read_events

from pylot.core.util.dataprocessing import Metadata
from pylot.core.util.obspyDMT_interface import qml_from_obspyDMT


def get_event_obspy_dmt(eventdir):
    event_pkl_file = os.path.join(eventdir, 'info', 'event.pkl')
    if not os.path.exists(event_pkl_file):
        raise IOError('Could not find event path for event: {}'.format(eventdir))
    event = qml_from_obspyDMT(event_pkl_file)
    return event


def get_event_pylot(eventdir, extension=''):
    event_id = get_event_id(eventdir)
    filename = os.path.join(eventdir, 'PyLoT_{}{}.xml'.format(event_id, extension))
    if not os.path.isfile(filename):
        return
    cat = read_events(filename)
    return cat[0]


def get_event_id(eventdir):
    event_id = os.path.split(eventdir)[-1]
    return event_id


def get_picks(eventdir, extension=''):
    event_id = get_event_id(eventdir)
    filename = 'PyLoT_{}{}.xml'
    filename = filename.format(event_id, extension)
    fpath = os.path.join(eventdir, filename)
    fpaths = glob.glob(fpath)
    if len(fpaths) == 1:
        cat = read_events(fpaths[0])
        picks = cat[0].picks
        return picks
    elif len(fpaths) == 0:
        print('get_picks: File not found: {}'.format(fpath))
        return
    print(f'WARNING: Ambiguous pick file specification. Found the following pick files {fpaths}\nFilemask: {fpath}')
    return


def write_json(object, fname):
    with open(fname, 'w') as outfile:
        json.dump(object, outfile, sort_keys=True, indent=4)


def get_metadata(eventdir):
    metadata_path = os.path.join(eventdir, 'resp')
    metadata = Metadata(inventory=metadata_path, verbosity=0)
    return metadata
