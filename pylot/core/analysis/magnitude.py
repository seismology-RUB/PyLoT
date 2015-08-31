# -*- coding: utf-8 -*-
"""
Created August/September 2015.

:author: Ludger Küperkoch / MAGS2 EP3 working group
"""

import matplotlib.pyplot as plt
import numpy as np
from obspy.core import Stream
from pylot.core.pick.utils import getsignalwin

class Magnitude(object):
    '''
    Superclass for calculating Wood-Anderson peak-to-peak
    amplitudes, local magnitudes and moment magnitudes.
    '''

    def __init__(self, wfstream, To, pwin, iplot):
        '''
        :param: wfstream
        :type: `~obspy.core.stream.Stream

        :param: To, onset time, P- or S phase
        :type:  float

        :param: pwin, pick window [To To+pwin] to get maximum
                peak-to-peak amplitude
        :type:  float

        :param: iplot, no. of figure window for plotting interims results
        :type: integer

        '''
        
        assert isinstance(wfstream, Stream), "%s is not a stream object" % str(wfstream)

        self.setwfstream(wfstream)
        self.setTo(To)
        self.setpwin(pwin)
        self.setiplot(iplot)
        self.calcwapp()


    def getwfstream(self):
        return self.wfstream

    def setwfstream(self, wfstream):
        self.wfstream = wfstream

    def getTo(self):
        return self.To

    def setTo(self, To):
        self.To = To

    def getpwin(self):
        return self.pwin

    def setpwin(self, pwin):
        self.pwin = pwin
   
    def getiplot(self):
        return self.iplot

    def setiplot(self, iplot):
        self.iplot = iplot

    def getwapp(self):
        return self.wapp

    def calcwapp(self):
        self.wapp = None

class WApp(Magnitude):
    '''
    Method to derive peak-to-peak amplitude as seen on a Wood-Anderson-
    seismograph. Has to be derived from corrected traces!
    '''

    def calcwapp(self):
        print "Getting Wood-Anderson peak-to-peak amplitude ..."
        print "Simulating Wood-Anderson seismograph ..."

        self.wapp = None
        stream = self.getwfstream()

        # poles, zeros and sensitivity of WA seismograph
        # (see Uhrhammer & Collins, 1990, BSSA, pp. 702-716)
        paz_wa = {
            'poles': [5.6089 - 5.4978j, -5.6089 - 5.4978j],
            'zeros': [0j, 0j],
            'gain': 2080,
            'sensitivity': 1} 

        stream.simulate(paz_remove=None, paz_simulate=paz_wa)

        trH1 = stream[0].data
        trH2 = stream[1].data
        ilen = min([len(trH1), len(trH2)])
        # get RMS of both horizontal components
        sqH = np.sqrt(np.power(trH1[0:ilen], 2) + np.power(trH2[0:ilen], 2))
        # get time array
        th = np.arange(0, len(sqH) * stream[0].stats.delta, stream[0].stats.delta)
        # get maximum peak within pick window
        iwin = getsignalwin(th, self.getTo(), self.getpwin())
        self.wapp = np.max(sqH[iwin])
        print "Determined Wood-Anderson peak-to-peak amplitude: %f mm" % self.wapp
        if self.getiplot() > 1:
            stream.plot()
            f = plt.figure(2)
            plt.plot(th, sqH)
            plt.plot(th[iwin], sqH[iwin], 'g')
            plt.plot([self.getTo(), self.getTo()], [0, max(sqH)], 'r', linewidth=2)
            plt.title('Station %s, RMS Horizontal Traces, WA-peak-to-peak=%4.1f mm' \
                      % (stream[0].stats.station, self.wapp))
            plt.xlabel('Time [s]')
            plt.ylabel('Displacement [mm]')
            plt.show()
            raw_input()
            plt.close(f)

    

