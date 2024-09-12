#!/bin/bash

#ulimit -s 8192
#ulimit -v $(ulimit -v | awk '{printf("%d",$1*0.95)}')
#ulimit -v

#655360

source /opt/anaconda3/etc/profile.d/conda.sh
conda activate pylot_311
NSLOTS=20

#qsub -l low -cwd -l "os=*stretch" -pe smp 40 submit_pick_corr_correction.sh
#$ -l low
#$ -l h_vmem=6G
#$ -cwd
#$ -pe smp 20
#$ -N corr_pick


export PYTHONPATH="$PYTHONPATH:/home/marcel/git/pylot_tools/"
export PYTHONPATH="$PYTHONPATH:/home/marcel/git/"
export PYTHONPATH="$PYTHONPATH:/home/marcel/git/pylot/"

#export MKL_NUM_THREADS=${NSLOTS:=1}
#export NUMEXPR_NUM_THREADS=${NSLOTS:=1}
#export OMP_NUM_THREADS=${NSLOTS:=1}

#python pick_correlation_correction.py '/data/AlpArray_Data/dmt_database_mantle_M5.8-6.0' '/home/marcel/.pylot/pylot_alparray_mantle_corr_stack_0.03-0.5.in' -pd -n ${NSLOTS:=1} -istart 0 -istop 100
#python pick_correlation_correction.py '/data/AlpArray_Data/dmt_database_mantle_M5.8-6.0' '/home/marcel/.pylot/pylot_alparray_mantle_corr_stack_0.03-0.5.in' -pd -n ${NSLOTS:=1} -istart 100 -istop 200
#python pick_correlation_correction.py '/data/AlpArray_Data/dmt_database_mantle_M6.0-6.5' '/home/marcel/.pylot/pylot_alparray_mantle_corr_stack_0.03-0.5.in' -pd -n ${NSLOTS:=1} -istart 0 -istop 100
#python pick_correlation_correction.py '/data/AlpArray_Data/dmt_database_mantle_M5.8-6.0' '/home/marcel/.pylot/pylot_alparray_mantle_corr_stack_0.03-0.5.in' -pd -n ${NSLOTS:=1} -istart 100 -istop 200
#python pick_correlation_correction.py 'H:\sciebo\dmt_database' 'H:\Sciebo\dmt_database\pylot_alparray_mantle_corr_S_0.01-0.2.in' -pd -n 4 -t

pylot_infile='/home/marcel/.pylot/pylot_alparray_syn_fwi_mk6_it3.in'
#pylot_infile='/home/marcel/.pylot/pylot_adriaarray_corr_P_and_S.in'

# THIS SCRIPT SHOLD BE CALLED BY "submit_to_grid_engine.py" using the following line:
python pick_correlation_correction.py $1 $pylot_infile -pd -n ${NSLOTS:=1} -istart $2 --params 'parameters_fwi_mk6_it3.yaml'
#--event_blacklist eventlist.txt
