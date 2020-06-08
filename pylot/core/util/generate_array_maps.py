#!/usr/bin/env python
# -*- coding: utf-8 -*-
# small script that creates array maps for each event within a previously generated PyLoT project

import os
import multiprocessing

from PyLoT import Project
from pylot.core.util.dataprocessing import Metadata
from pylot.core.util.array_map import Array_map

import matplotlib.pyplot as plt

def main(project_file_path, manual=False, auto=True, file_format='png', f_ext='', ncores=None):
    project = Project.load(project_file_path)
    nEvents = len(project.eventlist)
    input_list = []

    for index, event in enumerate(project.eventlist):
        # MP MP TESTING +++
        #if not eventdir.endswith('20170908_044946.a'):
        #    continue
        # MP MP ----
        kwargs = dict(project=project, event=event, nEvents=nEvents, index=index, manual=manual, auto=auto,
                      file_format=file_format, f_ext=f_ext)
        input_list.append(kwargs)

    if ncores == 1:
        for item in input_list:
            array_map_worker(item)
    else:
        pool = multiprocessing.Pool(ncores)
        result = pool.map(array_map_worker, input_list)
        pool.close()
        pool.join()

def array_map_worker(input_dict):
    event = input_dict['event']
    eventdir = event.path
    print('Working on event: {} ({}/{})'.format(eventdir, input_dict['index'] + 1, input_dict['nEvents']))
    # check for picks
    manualpicks = event.getPicks()
    autopicks = event.getAutopicks()
    # prepare event and get metadata
    metadata_path = os.path.join(eventdir, 'resp')
    metadata = None
    for pick_type in ['manual', 'auto']:
        if pick_type == 'manual' and (not manualpicks or not input_dict['manual']):
            continue
        if pick_type == 'auto' and (not autopicks or not input_dict['auto']):
            continue
        if not metadata:
            metadata = Metadata(inventory=metadata_path, verbosity=0)
        # create figure to plot on
        fig = plt.figure(figsize=(16,9))
        # create array map object
        map = Array_map(None, metadata, parameter=input_dict['project'].parameter, figure=fig,
                        width=2.13e6, height=1.2e6, pointsize=15., linewidth=1.0)
        # set combobox to auto/manual to plot correct pick type
        map.comboBox_am.setCurrentIndex(map.comboBox_am.findText(pick_type))
        # add picks to map and save file
        map.refresh_drawings(manualpicks, autopicks)
        fpath_out = os.path.join(eventdir, 'array_map_{}_{}{}.{}'.format(event.pylot_id, pick_type, input_dict['f_ext'],
                                                                         input_dict['file_format']))
        fig.savefig(fpath_out, dpi=300.)
        print('Wrote file: {}'.format(fpath_out))

if __name__ == '__main__':
    dataroot = '/home/marcel'
    infiles=['alparray_all_events_0.03-0.1_mantle_correlated_v3.plp']

    for infile in infiles:
        main(os.path.join(dataroot, infile), f_ext='_correlated_0.1Hz', ncores=10)
    #main('E:\Shared\AlpArray\\test_aa.plp', f_ext='_correlated_0.5Hz', ncores=1)
    #main('/home/marcel/alparray_m6.5-6.9_mantle_correlated_v3.plp', f_ext='_correlated_0.5Hz')
