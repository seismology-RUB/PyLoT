#!/bin/bash

#$ -l low
#$ -cwd
#$ -pe smp 40
#$ -l mem=2G
#$ -l h_vmem=2G
#$ -l os=*stretch

python ./autoPyLoT.py -i /home/marcel/.pylot/pylot_alparray_mantle_corr_stack_0.03-0.5.in -dmt processed -c $NSLOTS
