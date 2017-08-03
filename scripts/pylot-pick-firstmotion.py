#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
   Created Mar 2015
   Function to derive first motion (polarity) for given phase onset based on zero crossings.

   :author: MAGS2 EP3 working group / Ludger Kueperkoch
"""

import argparse

import obspy
from pylot.core.pick.utils import fmpicker

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--Xraw', type=obspy.core.stream.Stream,
                        help='unfiltered time series (seismogram) read with obspy module read')
    parser.add_argument('--Xfilt', type=obspy.core.stream.Stream,
                        help='filtered time series (seismogram) read with obspy module read')
    parser.add_argument('--pickwin', type=float, help='length of pick window [s] for first motion determination')
    parser.add_argument('--Pick', type=float, help='Onset time of most likely pick')
    parser.add_argument('--iplot', type=int, help='if set, figure no. iplot occurs')
    args = parser.parse_args()
    fmpicker(args.Xraw, args.Xfilt, args.pickwin, args.Pick, args.iplot)
