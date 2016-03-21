#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
   Created Mar 2015
   Transcription of the rezipe of Diehl et al. (2009) for consistent phase
   picking. For a given inital (the most likely) pick, the corresponding earliest
   and latest possible pick is calculated based on noise measurements in front of
   the most likely pick and signal wavelength derived from zero crossings.

   :author: Ludger Kueperkoch / MAGS2 EP3 working group
"""

import argparse
import obspy
from pylot.core.pick.utils import earllatepicker


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--X', type=~obspy.core.stream.Stream, help='time series (seismogram) read with obspy module read')
    parser.add_argument('--nfac', type=int, help='(noise factor), nfac times noise level to calculate latest possible pick')
    parser.add_argument('--TSNR', type=tuple, help='length of time windows around pick used to determine SNR \
                        [s] (Tnoise, Tgap, Tsignal)')
    parser.add_argument('--Pick1', type=float, help='Onset time of most likely pick')
    parser.add_argument('--iplot', type=int, help='if set, figure no. iplot occurs')
    args = parser.parse_args()
    earllatepicker(args.X, args.nfac, args.TSNR, args.Pick1, args.iplot)
