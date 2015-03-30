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
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import Stream
import argparse

def earllatepicker(X, nfac, TSNR, Pick1, iplot=None):
    '''
    Function to derive earliest and latest possible pick after Diehl & Kissling (2009) 
    as reasonable uncertainties. Latest possible pick is based on noise level,
    earliest possible pick is half a signal wavelength in front of most likely 
    pick given by PragPicker or manually set by analyst. Most likely pick
    (initial pick Pick1) must be given. 
    
    :param: X, time series (seismogram)
    :type:  `~obspy.core.stream.Stream`
   
    :param: nfac (noise factor), nfac times noise level to calculate latest possible pick
    :type: int

    :param: TSNR, length of time windows around pick used to determine SNR [s]
    :type: tuple (T_noise, T_gap, T_signal)
 
    :param: Pick1, initial (most likely) onset time, starting point for earllatepicker
    :type: float

    :param: iplot, if given, results are plotted in figure(iplot)
    :type: int
    '''

    assert isinstance(X, Stream), "%s is not a stream object" % str(X)

    LPick = None 
    EPick = None
    PickError = None
    if Pick1 is not None:
       print 'earllatepicker: Get earliest and latest possible pick relative to most likely pick ...'

       x =X[0].data
       t = np.arange(0, X[0].stats.npts / X[0].stats.sampling_rate, X[0].stats.delta)
       #some parameters needed:
       tnoise = TSNR[0]  #noise window length for calculating noise level
       tsignal = TSNR[2] #signal window length
       tsafety = TSNR[1] #safety gap between signal onset and noise window

       #get latest possible pick
       #get noise window
       inoise = np.where((t <= max([Pick1 - tsafety, 0])) \
                & (t >= max([Pick1 - tnoise - tsafety, 0])))    
       #get signal window
       isignal = np.where((t <= min([Pick1 + tsignal + tsafety, len(x)])) \
                 & (t >= Pick1))
       #calculate noise level
       nlevel = max(abs(x[inoise])) * nfac
       #get time where signal exceeds nlevel
       ilup = np.where(x[isignal] > nlevel)
       ildown = np.where(x[isignal] < -nlevel)
       if len(ilup[0]) <= 1 and len(ildown[0]) <= 1:
          print 'earllatepicker: Signal lower than noise level, misspick?'
          return
       il = min([ilup[0][0], ildown[0][0]])
       LPick = t[isignal][il]

       #get earliest possible pick
       #get next 2 zero crossings after most likely pick
       #initial onset is assumed to be the first zero crossing
       zc = []
       zc.append(Pick1)
       i = 0
       for j in range(isignal[0][1],isignal[0][len(t[isignal]) - 1]):
           i = i+ 1
           if x[j-1] <= 0 and x[j] >= 0:
              zc.append(t[isignal][i])
           elif x[j-1] > 0 and x[j] <= 0:
                zc.append(t[isignal][i])
           if len(zc) == 3:
               break
       #calculate maximum period of signal out of zero crossings
       Ts = max(np.diff(zc))
       #Ts/4 is assumed as time difference between most likely and earliest possible pick!
       EPick = Pick1 - Ts/4

       #get symmetric pick error as mean from earliest and latest possible pick
       #by weighting latest possible pick two times earliest possible pick
       diffti_tl = LPick -Pick1 
       diffti_te = Pick1 - EPick
       PickError = (diffti_te + 2 * diffti_tl)  / 3

       if iplot is not None:
          plt.figure(iplot)
          p1, = plt.plot(t, x, 'k')
          p2, = plt.plot(t[inoise], x[inoise])
          p3, = plt.plot(t[isignal], x[isignal], 'r')
          p4, = plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
          p5, = plt.plot(zc, [0, 0, 0], '*g', markersize=14)
          plt.legend([p1, p2, p3, p4, p5], ['Data', 'Noise Window', 'Signal Window', 'Noise Level', 'Zero Crossings'], \
                      loc='best')
          plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
          plt.plot([Pick1, Pick1], [max(x), -max(x)], 'b', linewidth=2)
          plt.plot([LPick, LPick], [max(x)/2, -max(x)/2], '--k')
          plt.plot([EPick, EPick], [max(x)/2, -max(x)/2], '--k')
          plt.plot([Pick1 + PickError, Pick1 + PickError], [max(x)/2, -max(x)/2], 'r--')
          plt.plot([Pick1 - PickError, Pick1 - PickError], [max(x)/2, -max(x)/2], 'r--')
          plt.xlabel('Time [s] since %s' % X[0].stats.starttime)
          plt.yticks([])
          ax = plt.gca()
          ax.set_xlim([t[inoise[0][0]] - 2, t[isignal[0][len(isignal) - 1]] + 3])
          plt.title('Earliest-/Latest Possible/Most Likely Pick & Symmetric Pick Error, %s' % X[0].stats.station)
          plt.show()
          raw_input()
          plt.close(iplot)

    elif Pick1 == None: 
         print 'earllatepicker: No initial onset time given! Check input!'
         return

    return EPick, LPick, PickError

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
