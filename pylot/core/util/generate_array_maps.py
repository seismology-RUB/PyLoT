#!/usr/bin/env python
# -*- coding: utf-8 -*-
# small script that creates array maps for each event within a previously generated PyLoT project

import os
import multiprocessing
import sys
import glob
import matplotlib
matplotlib.use('Qt5Agg')
sys.path.append(os.path.join('/'.join(sys.argv[0].split('/')[:-1]), '../../..'))

from PyLoT import Project
from pylot.core.util.dataprocessing import Metadata
from pylot.core.util.array_map import Array_map

import matplotlib.pyplot as plt
import argparse


def main(project_file_path, manual=False, auto=True, file_format='png', f_ext='', ncores=None):
    project = Project.load(project_file_path)
    nEvents = len(project.eventlist)
    input_list = []
    print('\n')
    for index, event in enumerate(project.eventlist):
        kwargs = dict(project=project, event=event, nEvents=nEvents, index=index, manual=manual, auto=auto,
                      file_format=file_format, f_ext=f_ext)
        input_list.append(kwargs)

    if ncores == 1:
        for item in input_list:
            array_map_worker(item)
    else:
        pool = multiprocessing.Pool(ncores)
        pool.map(array_map_worker, input_list)
        pool.close()
        pool.join()


def array_map_worker(input_dict):
    event = input_dict['event']
    eventdir = event.path
    print('Working on event: {} ({}/{})'.format(eventdir, input_dict['index'] + 1, input_dict['nEvents']))
    xml_picks = glob.glob(os.path.join(eventdir, f'*{input_dict["f_ext"]}.xml'))
    if not len(xml_picks):
        print('Event {} does not have any picks associated with event file extension {}'. format(eventdir, input_dict['f_ext']))
        return
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
        fig, ax = plt.subplots(figsize=(16, 9))
        # create array map object
        map = Array_map(None, metadata, parameter=input_dict['project'].parameter, axes=ax,
                        width=2.13e6, height=1.2e6, pointsize=25., linewidth=1.0)
        # set combobox to auto/manual to plot correct pick type
        map.comboBox_am.setCurrentIndex(map.comboBox_am.findText(pick_type))
        # add picks to map and save file
        map.refresh_drawings(manualpicks, autopicks)
        fpath_out = os.path.join(eventdir, 'array_map_{}_{}{}.{}'.format(event.pylot_id, pick_type, input_dict['f_ext'],
                                                                         input_dict['file_format']))
        fig.savefig(fpath_out, dpi=300.)
        print('Wrote file: {}'.format(fpath_out))


if __name__ == '__main__':
    cl = argparse.ArgumentParser()
    cl.add_argument('--dataroot', help='Directory containing the PyLoT .plp file', type=str)
    cl.add_argument('--infiles', help='.plp files to use', nargs='+')
    cl.add_argument('--ncores', hepl='Specify number of parallel processes', type=int, default=1)
    args = cl.parse_args()

    for infile in args.infiles:
        main(os.path.join(args.dataroot, infile), f_ext='_correlated_0.03-0.1', ncores=args.ncores)

