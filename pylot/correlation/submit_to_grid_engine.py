#!/usr/bin/env python

import subprocess

fnames = [
    ('/data/AlpArray_Data/dmt_database_synth_model_mk6_it3_no_rotation', 0),
          ]

#fnames = [('/data/AlpArray_Data/dmt_database_mantle_0.01-0.2_SKS-phase', 0),
#          ('/data/AlpArray_Data/dmt_database_mantle_0.01-0.2_S-phase', 0),]

####
script_location = '/home/marcel/VersionCtrl/git/code_base/correlation_picker/submit_pick_corr_correction.sh'
####

for fnin, istart in fnames:
    input_cmds = f'qsub -q low.q@minos15,low.q@minos14,low.q@minos13,low.q@minos12,low.q@minos11 {script_location} {fnin} {istart}'

    print(input_cmds)
    print(subprocess.check_output(input_cmds.split()))



