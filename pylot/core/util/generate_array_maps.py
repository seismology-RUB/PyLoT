#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from PyLoT import Project
from pylot.core.util.dataprocessing import Metadata
from pylot.core.util.array_map import Array_map

import matplotlib.pyplot as plt

def main(project_file_path, manual=False, auto=True):
    project = Project.load(project_file_path)
    nEvents = len(project.eventlist)

    for index, event in enumerate(project.eventlist):
        eventdir = event.path
        print('Working on event: {} ({}/{})'.format(eventdir, index + 1, nEvents))
        # check for picks
        manualpicks = event.getPicks()
        autopicks = event.getAutopicks()
        # prepare event and get metadata
        metadata_path = os.path.join(eventdir, 'resp')
        metadata = Metadata(inventory=metadata_path, verbosity=0)
        for pick_type in ['manual', 'auto']:
            if pick_type == 'manual' and (not manualpicks or not manual):
                continue
            if pick_type == 'auto' and (not autopicks or not auto):
                continue
            # create figure to plot on
            fig = plt.figure(figsize=(16,9))
            # create array map object
            map = Array_map(None, metadata, figure=fig, width=2.13e6, height=1.2e6)
            # set combobox to auto/manual to plot correct pick type
            map.comboBox_am.setCurrentIndex(map.comboBox_am.findText(pick_type))
            # add picks to map and save file
            map.refresh_drawings(manualpicks, autopicks)
            fpath_out = os.path.join(eventdir, 'array_map_{}_{}.svg'.format(event.pylot_id, pick_type))
            fig.savefig(fpath_out, dpi=300.)
            print('Wrote file: {}'.format(fpath_out))

if __name__ == '__main__':
    main('/home/marcel/pylot_m7_mantle_correlated.plp')
