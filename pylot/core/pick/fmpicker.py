#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
   Created Mar 2015
   Function to derive first motion (polarity) for given phase onset based on zero crossings.

   :author: MAGS2 EP3 working group / Ludger Kueperkoch	
"""
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import Stream
import argparse

def fmpicker(Xraw, Xfilt, pickwin, Pick, iplot=None):
    '''
    Function to derive first motion (polarity) of given phase onset Pick.
    Calculation is based on zero crossings determined within time window pickwin
    after given onset time. 
    
    :param: Xraw, unfiltered time series (seismogram)
    :type:  `~obspy.core.stream.Stream`
   
    :param: Xfilt, filtered time series (seismogram)
    :type:  `~obspy.core.stream.Stream`

    :param: pickwin, time window after onset Pick within zero crossings are calculated
    :type: float

    :param: Pick, initial (most likely) onset time, starting point for fmpicker
    :type: float

    :param: iplot, if given, results are plotted in figure(iplot)
    :type: int
    '''

    assert isinstance(Xraw, Stream), "%s is not a stream object" % str(Xraw)
    assert isinstance(Xfilt, Stream), "%s is not a stream object" % str(Xfilt)

    FM = None
    if Pick is not None:
       print 'fmpicker: Get first motion (polarity) of onset using unfiltered seismogram...'

       xraw = Xraw[0].data
       xfilt = Xfilt[0].data
       t = np.arange(0, Xraw[0].stats.npts / Xraw[0].stats.sampling_rate, Xraw[0].stats.delta)
       #get pick window
       ipick = np.where((t <= min([Pick + pickwin, len(Xraw[0])])) & (t >= Pick))
       #remove mean
       xraw[ipick] = xraw[ipick] - np.mean(xraw[ipick])
       xfilt[ipick] = xfilt[ipick] - np.mean(xfilt[ipick])

       #get next zero crossing after most likely pick
       #initial onset is assumed to be the first zero crossing
       #first from unfiltered trace
       zc1 = []
       zc1.append(Pick)
       index1 = []
       i = 0
       for j in range(ipick[0][1],ipick[0][len(t[ipick]) - 1]):
           i = i+ 1
           if xraw[j-1] <= 0 and xraw[j] >= 0:
              zc1.append(t[ipick][i])
              index1.append(i)
           elif xraw[j-1] > 0 and xraw[j] <= 0:
              zc1.append(t[ipick][i])
              index1.append(i)
           if len(zc1) == 3:
               break

       #if time difference betweeen 1st and 2cnd zero crossing
       #is too short, get time difference between 1st and 3rd
       #to derive maximum 
       if zc1[1] - zc1[0] <= Xraw[0].stats.delta:
          li1 = index1[1]
       else:
          li1 = index1[0]
       if np.size(xraw[ipick[0][1]:ipick[0][li1]]) == 0:
          print 'earllatepicker: Onset on unfiltered trace too emergent for first motion determination!'
          P1 = None
       else: 
          imax1 = np.argmax(abs(xraw[ipick[0][1]:ipick[0][li1]]))
          islope1 = np.where((t >= Pick) & (t <= Pick + t[imax1]))
          #calculate slope as polynomal fit of order 1
          xslope1 = np.arange(0, len(xraw[islope1]), 1)
          P1 = np.polyfit(xslope1, xraw[islope1], 1)
          datafit1 = np.polyval(P1, xslope1)

       #now using filterd trace
       #next zero crossing after most likely pick
       zc2 = []
       zc2.append(Pick)
       index2 = []
       i = 0
       for j in range(ipick[0][1],ipick[0][len(t[ipick]) - 1]):
           i = i+ 1
           if xfilt[j-1] <= 0 and xfilt[j] >= 0:
              zc2.append(t[ipick][i])
              index2.append(i)
           elif xfilt[j-1] > 0 and xfilt[j] <= 0:
              zc2.append(t[ipick][i])
              index2.append(i)
           if len(zc2) == 3:
               break

       #if time difference betweeen 1st and 2cnd zero crossing
       #is too short, get time difference between 1st and 3rd
       #to derive maximum 
       if zc2[1] - zc2[0] <= Xfilt[0].stats.delta:
          li2 = index2[1]
       else:
          li2 = index2[0]
       if np.size(xfilt[ipick[0][1]:ipick[0][li2]]) == 0:
          print 'earllatepicker: Onset on filtered trace too emergent for first motion determination!'
          P2 = None
       else: 
          imax2 = np.argmax(abs(xfilt[ipick[0][1]:ipick[0][li2]]))
          islope2 = np.where((t >= Pick) & (t <= Pick + t[imax2]))
          #calculate slope as polynomal fit of order 1
          xslope2 = np.arange(0, len(xfilt[islope2]), 1)
          P2 = np.polyfit(xslope2, xfilt[islope2], 1)
          datafit2 = np.polyval(P2, xslope2)
       
       #compare results
       if P1 is not None and P2 is not None:
          if P1[0] < 0 and P2[0] < 0:
             FM = 'D'
          elif P1[0] >= 0 and P2[0] < 0:
             FM = '-'
          elif P1[0] < 0 and P2[0]>= 0:
             FM = '-'
          elif P1[0] > 0 and P2[0] > 0:
             FM = 'U'
          elif P1[0] <= 0 and P2[0] > 0:
             FM = '+'
          elif P1[0] > 0 and P2[0] <= 0:
             FM = '+'

    if iplot is not None:
       plt.figure(iplot)
       plt.subplot(2,1,1)
       plt.plot(t, xraw, 'k')
       p1, = plt.plot([Pick, Pick], [max(xraw), -max(xraw)], 'b', linewidth=2)
       if P1 is not None:
          p2, = plt.plot(t[islope1], xraw[islope1])
          p3, = plt.plot(zc1, np.zeros(len(zc1)), '*g', markersize=14)
          p4, = plt.plot(t[islope1], datafit1, '--g', linewidth=2) 
          plt.legend([p1, p2, p3, p4], ['Pick', 'Slope Window', 'Zero Crossings', 'Slope'], \
                      loc='best')
          plt.text(Pick + 0.02, max(xraw) / 2, '%s' % FM, fontsize=14) 
          ax = plt.gca()
          ax.set_xlim([t[islope1[0][0]] - 0.1, t[islope1[0][len(islope1) - 1]] + 0.3])
       plt.yticks([])
       plt.title('First-Motion Determination, %s, Unfiltered Data' % Xraw[0].stats.station)

       plt.subplot(2,1,2)
       plt.title('First-Motion Determination, Filtered Data')
       plt.plot(t, xfilt, 'k')
       p1, = plt.plot([Pick, Pick], [max(xfilt), -max(xfilt)], 'b', linewidth=2)
       if P2 is not None:
          p2, = plt.plot(t[islope2], xfilt[islope2])
          p3, = plt.plot(zc2, np.zeros(len(zc2)), '*g', markersize=14)
          p4, = plt.plot(t[islope2], datafit2, '--g', linewidth=2) 
          plt.text(Pick + 0.02, max(xraw) / 2, '%s' % FM, fontsize=14) 
          ax = plt.gca()
          ax.set_xlim([t[islope2[0][0]] - 0.1, t[islope2[0][len(islope2) - 1]] + 0.3])
       plt.xlabel('Time [s] since %s' % Xraw[0].stats.starttime)
       plt.yticks([])
       plt.show()
       raw_input()
       plt.close(iplot)

    return FM

if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('--Xraw', type=~obspy.core.stream.Stream, help='unfiltered time series (seismogram) read with obspy module read')
   parser.add_argument('--Xfilt', type=~obspy.core.stream.Stream, help='filtered time series (seismogram) read with obspy module read')
   parser.add_argument('--pickwin', type=float, help='length of pick window [s] for first motion determination')
   parser.add_argument('--Pick', type=float, help='Onset time of most likely pick')
   parser.add_argument('--iplot', type=int, help='if set, figure no. iplot occurs')
   args = parser.parse_args()
   earllatepicker(args.Xraw, args.Xfilt, args.pickwin, args.Pick, args.iplot)

