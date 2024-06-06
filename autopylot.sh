#!/bin/bash

#$ -l low
#$ -cwd
#$ -pe smp 40
##$ -l mem=3G
#$ -l h_vmem=6G
#$ -l os=*stretch

conda activate pylot_311

python ./autoPyLoT.py -i /home/marcel/.pylot/pylot_adriaarray.in -c 20 -dmt processed
