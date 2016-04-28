#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
   Created Mar/Apr 2015
   Function to calculate SNR of certain part of seismogram relative
   to given time. Returns SNR and SNR [dB].

   :author: Ludger Kueperkoch /MAGS EP3 working group
"""

import argparse
import obspy
from pylot.core.pick.utils import getSNR

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', '-d', type=obspy.core.stream.Stream,
                        help='time series (seismogram) read with obspy module '
                             'read',
                        dest='data')
    parser.add_argument('--tsnr', '-s', type=tuple,
                        help='length of time windows around pick used to '
                             'determine SNR [s] (Tnoise, Tgap, Tsignal)',
                        dest='tsnr')
    parser.add_argument('--time', '-t', type=float,
                        help='initial time from which noise and signal windows '
                             'are calculated',
                        dest='time')
    args = parser.parse_args()
    print getSNR(args.data, args.tsnr, args.time)
