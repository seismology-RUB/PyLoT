#!/bin/bash

#$ -l low
#$ -cwd
#$ -pe smp 40
#$ -l mem=2G
#$ -l h_vmem=2G
#$ -l os=*stretch

conda activate pylot_38

python ./autoPyLoT.py -i /home/marcel/.pylot/pylot_adriaarray.in -c 20 -dmt processed
