#!/usr/bin/env python

import subprocess

fnames = [
    ('/data/AdriaArray_Data/dmt_database_mantle_M5.0-5.4', 0),
    ('/data/AdriaArray_Data/dmt_database_mantle_M5.4-5.7', 0),
    ('/data/AdriaArray_Data/dmt_database_mantle_M5.7-6.0', 0),
    ('/data/AdriaArray_Data/dmt_database_mantle_M6.0-6.3', 0),
    ('/data/AdriaArray_Data/dmt_database_mantle_M6.3-10.0', 0),
    # ('/data/AdriaArray_Data/dmt_database_ISC_mantle_M5.0-5.4', 0),
    # ('/data/AdriaArray_Data/dmt_database_ISC_mantle_M5.4-5.7', 0),
    # ('/data/AdriaArray_Data/dmt_database_ISC_mantle_M5.7-6.0', 0),
    # ('/data/AdriaArray_Data/dmt_database_ISC_mantle_M6.0-10.0', 0),
          ]

#fnames = [('/data/AlpArray_Data/dmt_database_mantle_0.01-0.2_SKS-phase', 0),
#          ('/data/AlpArray_Data/dmt_database_mantle_0.01-0.2_S-phase', 0),]

####
script_location = '/home/marcel/VersionCtrl/git/pylot/pylot/correlation/submit_pick_corr_correction.sh'
####

for fnin, istart in fnames:
    input_cmds = f'qsub -q low.q@minos15,low.q@minos14,low.q@minos13,low.q@minos12,low.q@minos11 {script_location} {fnin} {istart}'

    print(input_cmds)
    print(subprocess.check_output(input_cmds.split()))