#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
   Created Mar/Apr 2015
   Function to calculate SNR of certain part of seismogram relativ 
   to given time. Returns SNR and SNR [dB].

   :author: Ludger Kueperkoch /MAGS EP3 working group
"""
from obspy.core import Stream
import numpy as np

def getSNR(X, TSNR, t1):
    '''
    Function to calculate SNR of certain part of seismogram relative to
    given time (onset) out of given noise and signal windows. A safety gap
    between noise and signal part can be set. Returns SNR and SNR [dB].

    :param: X, time series (seismogram)
    :type:  `~obspy.core.stream.Stream`

    :param: TSNR, length of time windows [s] around t1 (onset) used to determine SNR
    :type: tuple (T_noise, T_gap, T_signal)

    :param: t1, initial time (onset) from which noise and signal windows are calculated
    :type: float
    '''

    assert isinstance(X, Stream), "%s is not a stream object" % str(X)

    SNR = None
    SNRdB = None
    x = X[0].data
    t = np.arange(0, X[0].stats.npts / X[0].stats.sampling_rate, X[0].stats.delta)
    #some parameters needed:
    tnoise = TSNR[0]  #noise window length for calculating noise level
    tsignal = TSNR[2] #signal window length
    tsafety = TSNR[1] #safety gap between signal onset and noise window
    #get noise window
    inoise = np.where((t <= max([t1 - tsafety, 0])) \
             & (t >= max([t1 - tnoise - tsafety, 0])))
    #get signal window
    isignal = np.where((t <= min([t1 + tsignal + tsafety, len(x)])) \
              & (t >= t1))
    if np.size(inoise) < 1:
       print 'getSNR: Empty array inoise, check noise window!'
       return
    elif np.size(isignal) < 1:
       print 'getSNR: Empty array isignal, check signal window!'
       return

    #calculate ratios
    SNR = max(abs(x[isignal])) / np.mean(abs(x[inoise]))
    SNRdB = 20 * np.log10(SNR)
    
    return SNR, SNRdB

if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('--X', type=~obspy.core.stream.Stream, help='time series (seismogram) read with obspy module read')
   parser.add_argument('--TSNR', type=tuple, help='length of time windows around pick used to determine SNR \
                       [s] (Tnoise, Tgap, Tsignal)')
   parser.add_argument('--t1', type=float, help='initial time from which noise and signal windows are calculated')
   args = parser.parse_args()
   getSNR(args.X, args.TSNR, args.t1)

