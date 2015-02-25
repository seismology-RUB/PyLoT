# -*- coding: utf-8 -*-
"""
Created Dec 2014 to Feb 2015
Implementation of the automated picking algorithms published and described in:

Kueperkoch, L., Meier, T., Lee, J., Friederich, W., & Egelados Working Group, 2010:
Automated determination of P-phase arrival times at regional and local distances
using higher order statistics, Geophys. J. Int., 181, 1159-1170

Kueperkoch, L., Meier, T., Bruestle, A., Lee, J., Friederich, W., & Egelados
Working Group, 2012: Automated determination of S-phase arrival times using
autoregressive prediction: application ot local and regional distances, Geophys. J. Int.,
188, 687-702.

The picks with the above described algorithms are assumed to be the most likely picks.
For each most likely pick the corresponding earliest and latest possible picks are
calculated after Diehl & Kissling (2009).

:author: MAGS2 EP3 working group / Ludger Kueperkoch	
"""
import numpy as np
import matplotlib.pyplot as plt
from pylot.core.pick.CharFuns import CharacteristicFunction

class AutoPicking(object):
    '''
    Superclass of different, automated picking algorithms applied on a CF determined
    using AIC, HOS, or AR prediction.
    ''' 
    def __init__(self, cf, nfac, TSNR, PickWindow, iplot=None, aus=None, Tsmooth=None, Pick1=None):
        '''
        :param: cf, characteristic function, on which the picking algorithm is applied
        :type: `~pylot.core.pick.CharFuns.CharacteristicFunction` object

        :param: nfac (noise factor), nfac times noise level to calculate latest possible pick
         in EarlLatePick
        :type: int

        :param: TSNR, length of time windows around pick used to determine SNR [s]
        :type: tuple (T_noise, T_gap, T_signal)

        :param: PickWindow, length of pick window [s]
        :type: float

        :param: iplot, no. of figure window for plotting interims results
        :type: integer

        :param: aus ("artificial uplift of samples"), find local minimum at i if aic(i-1)*(1+aus) >= aic(i)
        :type: float

        :param: Tsmooth, length of moving smoothing window to calculate smoothed CF [s]
        :type: float

        :param: Pick1, initial (prelimenary) onset time, starting point for PragPicker and
         EarlLatePick
        :type: float

        '''

        assert isinstance(cf, CharacteristicFunction), "%s is not a CharacteristicFunction object" % str(cf)

        self.cf = cf.getCF()
        self.Tcf = cf.getTimeArray()
        self.Data = cf.getXCF()
        self.dt = cf.getIncrement()
        self.setnfac(nfac)
        self.setTSNR(TSNR)
        self.setPickWindow(PickWindow)
        self.setiplot(iplot)
        self.setaus(aus)
        self.setTsmooth(Tsmooth)
        self.setpick1(Pick1)
        self.calcPick()

    def __str__(self):
        return '''\n\t{name} object:\n
        nfac:\t{nfac}\n
        TSNR:\t\t\t{TSNR}\n
        PickWindow:\t{PickWindow}\n
        aus:\t{aus}\n
        Tsmooth:\t{Tsmooth}\n
        Pick1:\t{Pick1}\n
        '''.format(name=type(self).__name__,
                   nfac=self.getnfac(),
                   TSNR=self.getTSNR(),
                   PickWindow=self.getPickWindow(),
                   aus=self.getaus(),
                   Tsmooth=self.getTsmooth(),
                   Pick1=self.getpick1())     

    def getnfac(self):
        return self.nfac 
              
    def setnfac(self, nfac):
        self.nfac = nfac
 
    def getTSNR(self):
        return self.TSNR

    def setTSNR(self, TSNR):
        self.TSNR = TSNR

    def getPickWindow(self):
        return self.PickWindow

    def setPickWindow(self, PickWindow):
        self.PickWindow = PickWindow

    def getaus(self):
        return self.aus

    def setaus(self, aus):
        self.aus = aus

    def setTsmooth(self, Tsmooth):
        self.Tsmooth = Tsmooth

    def getTsmooth(self):
        return self.Tsmooth

    def getpick(self):
        return self.Pick

    def getLpick(self):
        return self.LPick

    def getEpick(self):
        return self.EPick

    def getiplot(self):
        return self.iplot

    def setiplot(self, iplot):
        self.iplot = iplot

    def getpick1(self):
        return self.Pick1

    def setpick1(self, Pick1):
        self.Pick1 = Pick1

    def calcPick(self):
        self.Pick = None


class AICPicker(AutoPicking):
    '''
    Method to derive onset time of arriving phase based on CF
    derived from AIC.
    '''

    def calcPick(self):
 
        print 'Get onset time (pick) from AIC-CF ...'

        self.Pick = None
        #find NaN's
        nn = np.isnan(self.cf)
        if len(nn) > 1:
           self.cf[nn] = 0
        #taper AIC-CF to get rid off side maxima
        tap = np.hanning(len(self.cf))
        aic = tap * self.cf + max(abs(self.cf))
        #smooth AIC-CF
        ismooth = int(round(self.Tsmooth / self.dt))
        aicsmooth = np.zeros(len(aic))
        if len(aic) < ismooth:
           print 'AICPicker: Tsmooth larger than CF!'
           return
        else:
           for i in range(1, len(aic)):
              if i > ismooth:
                 ii1 = i - ismooth;
                 aicsmooth[i] = aicsmooth[i - 1] + (aic[i] - aic[ii1]) / ismooth
              else:
                 aicsmooth[i] = np.mean(aic[1 : i])
        #remove offset
        offset = abs(min(aic) - min(aicsmooth))
        aicsmooth = aicsmooth - offset
        #get maximum of CF as starting point
        icfmax = np.argmax(aic)
        
        #find minimum in front of maximum
        lpickwindow = int(round(self.PickWindow / self.dt))
        for i in range(icfmax - 1, max([icfmax - lpickwindow, 2]), -1): 
           if aicsmooth[i - 1] >= aicsmooth[i]:
              self.Pick = self.Tcf[i]
              break
        
        if self.iplot is not None:
           plt.figure(self.iplot)
           x = self.Data[0].data
           p1, = plt.plot(self.Tcf, x / max(x), 'k')
           p2, = plt.plot(self.Tcf, aicsmooth / max(aicsmooth), 'r')
           p3, = plt.plot([self.Pick, self.Pick], [-1 , 1], 'b', linewidth=2)
           plt.legend([p1, p2, p3], ['(HOS-/AR-) Data', 'Smoothed AIC-CF', 'AIC-Pick'])
           plt.xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
           plt.yticks([])
           plt.title(self.Data[0].stats.station)
           plt.show()
           raw_input()
           plt.close(self.iplot)

        if self.Pick == None:
           print 'AICPicker: Could not find minimum, picking window too short?'
  
        return self.Pick


class PragPicker(AutoPicking):
    '''
    Method of pragmatic picking exploiting information given by CF.
    '''

    def calcPick(self):

        if self.getpick1() is not None:
           print 'Get most likely pick from HOS- or AR-CF using pragmatic picking algorithm ...'

           self.Pick = None 
           #smooth CF
           ismooth = int(round(self.Tsmooth / self.dt))
           cfsmooth = np.zeros(len(self.cf))
           if len(self.cf) < ismooth:
              print 'PragPicker: Tsmooth larger than CF!'
              return
           else:
              for i in range(1, len(self.cf)):
                 if i > ismooth:
                    ii1 = i - ismooth;
                    cfsmooth[i] = cfsmooth[i - 1] + (self.cf[i] - self.cf[ii1]) / ismooth
                 else:
                    cfsmooth[i] = np.mean(self.cf[1 : i])

           #select picking window
           #which is centered around tpick1
           ipick = np.where((self.Tcf >= self.getpick1() - self.PickWindow / 2) \
                 & (self.Tcf <= self.getpick1() + self.PickWindow / 2))
           cfipick = self.cf[ipick] - np.mean(self.cf[ipick])
           Tcfpick = self.Tcf[ipick]
           cfsmoothipick = cfsmooth[ipick]- np.mean(self.cf[ipick])
           ipick1 = np.argmin(abs(self.Tcf - self.getpick1()))
           cfpick1 = 2 * self.cf[ipick1]

           #check trend of CF, i.e. differences of CF and adjust aus regarding this trend
           #prominent trend: decrease aus
           #flat: use given aus
           cfdiff = np.diff(cfipick);
           i0diff = np.where(cfdiff > 0)
           cfdiff = cfdiff[i0diff]
           minaus = min(cfdiff * (1 + self.aus));
           aus1 = max([minaus, self.aus]);

           #at first we look to the right until the end of the pick window is reached
           flagpick_r = 0
           flagpick_l = 0
           flagpick = 0
           lpickwindow = int(round(self.PickWindow / self.dt))
           for i in range(max(np.insert(ipick, 0, 2)), min([ipick1 + lpickwindow + 1, len(self.cf) - 1])):
              if self.cf[i + 1] > self.cf[i] and self.cf[i - 1] >= self.cf[i]:
                 if cfsmooth[i - 1] * (1 + aus1) >= cfsmooth[i]:
                    if cfpick1 >= self.cf[i]:     
                       pick_r = self.Tcf[i]
                       self.Pick = pick_r
                       flagpick_l = 1
                       cfpick_r = self.cf[i]
                       break

           #now we look to the left
           for i in range(ipick1, max([ipick1 - lpickwindow + 1, 2]), -1):
              if self.cf[i + 1] > self.cf[i] and self.cf[i - 1] >= self.cf[i]:
                 if cfsmooth[i - 1] * (1 + aus1) >= cfsmooth[i]:
                    if cfpick1 >= self.cf[i]:     
                       pick_l = self.Tcf[i]
                       self.Pick = pick_l
                       flagpick_r = 1
                       cfpick_l = self.cf[i]
                       break

           #now decide which pick: left or right?
           if flagpick_l > 0 and flagpick_r > 0 and cfpick_l <= cfpick_r:
              self.Pick = pick_l
           elif flagpick_l > 0 and flagpick_r > 0 and cfpick_l >= cfpick_r:
              self.Pick = pick_r

           if self.getiplot() is not None:
              plt.figure(self.getiplot())
              p1, = plt.plot(Tcfpick,cfipick, 'k')
              p2, = plt.plot(Tcfpick,cfsmoothipick, 'r')
              p3, = plt.plot([self.Pick, self.Pick], [min(cfipick), max(cfipick)], 'b', linewidth=2)
              plt.legend([p1, p2, p3], ['CF', 'Smoothed CF', 'Pick']) 
              plt.xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
              plt.yticks([])
              plt.title(self.Data[0].stats.station)
              plt.show()
              raw_input()
              plt.close(self.getiplot())

        else: 
           self.Pick = None
           print 'PragPicker: No initial onset time given! Check input!'
           return


class EarlLatePicker(AutoPicking):
    '''
    Method to derive earliest and latest possible pick after Diehl & Kissling (2009) 
    as reasonable uncertainties. Latest possible pick is based on noise level,
    earliest possible pick is half a signal wavelength in front of most likely 
    pick given by PragPicker. Most likely pick (initial pick) must be given. 
    '''
       
    def calcPick(self):

        self.LPick = None 
        self.EPick = None
        if self.getpick1() is not None:
           print 'Get earliest and latest possible pick relative to most likely pick ...'

           ti = self.getpick1()
           x = self.Data
           t = self.Tcf
           #some parmaters needed:
           tnoise = self.TSNR[0]  #noise window length for calculating noise level
           tsignal = self.TSNR[2] #signal window length
           tsafety = self.TSNR[1] #safety gap between signal onset and noise window

           #get latest possible pick
           #get noise window
           inoise = np.where((self.Tcf <= ti - tsafety) & (self.Tcf >= ti - tnoise - tsafety))    
           #get signal window
           isignal = np.where((self.Tcf <= ti + tsignal) & (self.Tcf >= ti))
           #calculate noise level
           if len(x) == 1:
              nlevel = max(abs(x[0].data[inoise])) * self.nfac
              #get time where signal exceeds nlevel
              ilup = np.where(x[0].data[isignal] > nlevel)
              ildown = np.where(x[0].data[isignal] < -nlevel)
              if len(ilup[0]) <= 1 and len(ildown[0]) <= 1:
                 print 'EarlLatePick: Signal lower than noise level, misspick?'
                 return
              il = min([ilup[0][0], ildown[0][0]])
              self.LPick = t[isignal][il]
           elif len(x) == 2:
              nlevel = max(np.sqrt(np.power(x[0].data[inoise], 2) + np.power(x[1].data[inoise], 2)))
              #get earliest time where signal exceeds nlevel
              ilup1 = np.where(x[0].data[isignal] > nlevel)
              ilup2 = np.where(x[1].data[isignal] > nlevel)
              ildown1 = np.where(x[0].data[isignal] < -nlevel)
              ildown2 = np.where(x[1].data[isignal] < -nlevel)
              ilup = min([ilup1[0][0], ilup2[0][0]])
              ildown = min([ildown1[0][0], ildown2[0][0]])
              if np.size(ilup) < 1 and np.size(ildown) < 1:
                 print 'EarlLatePick: Signal lower than noise level, misspick?'
                 return
              il = min([ilup, ildown])
              self.LPick = t[isignal][il]
           elif len(x) == 3:
              nlevel = max(np.sqrt(np.power(x[0].data[inoise], 2) + np.power(x[1].data[inoise], 2) + \
                       np.power(x[2].data[inoise], 2)))
              #get earliest time where signal exceeds nlevel
              ilup1 = np.where(x[0].data[isignal] > nlevel)
              ilup2 = np.where(x[1].data[isignal] > nlevel)
              ilup3 = np.where(x[2].data[isignal] > nlevel)
              ildown1 = np.where(x[0].data[isignal] < -nlevel)
              ildown2 = np.where(x[1].data[isignal] < -nlevel)
              ildown3 = np.where(x[2].data[isignal] < -nlevel)
              ilup = min([ilup1[0][0], ilup2[0][0], ilup3[0][0]])
              ildown = min([ildown1[0][0], ildown2[0][0], ildown3[0][0]])
              if np.size(ilup) < 1 and np.size(ildown) < 1:
                 print 'EarlLatePick: Signal lower than noise level, misspick?'
                 return
              il = min([ilup, ildown])
              self.LPick = t[isignal][il]

           #get earliest possible pick
           #get next 2 zero crossings after most likely pick
           #if there is one trace in stream
           if len(x) == 1:
              zc = []
              zc.append(ti)
              i = 0
              for j in range(isignal[0][1],isignal[0][len(t[isignal]) - 1]):
                  i = i+ 1
                  if x[0].data[j-1] <= 0 and x[0].data[j] >= 0:
                     zc.append(t[isignal][i])
                  elif x[0].data[j-1] > 0 and x[0].data[j] <= 0:
                     zc.append(t[isignal][i])
                  if len(zc) == 3:
                     break
              #calculate maximum period of signal out of zero crossings
              Ts = max(np.diff(zc))
           #if there are two traces in stream
           #get maximum of two signal periods
           if len(x) == 2:
              zc1 = []
              zc2 = []
              zc1.append(ti)
              zc2.append(ti)
              i = 0
              for j in range(isignal[0][1],isignal[0][len(t[isignal]) - 1]):
                  i = i+ 1
                  if x[0].data[j-1] <= 0 and x[0].data[j] >= 0:
                     zc1.append(t[isignal][i])
                  elif x[0].data[j-1] > 0 and x[0].data[j] <= 0:
                     zc1.append(t[isignal][i])
                  if x[1].data[j-1] <= 0 and x[1].data[j] >= 0:
                     zc2.append(t[isignal][i])
                  elif x[1].data[j-1] > 0 and x[1].data[j] <= 0:
                     zc2.append(t[isignal][i])
                  if len(zc1) >= 3 and len(zc2) >= 3:
                     break
              Ts = max([max(np.diff(zc1)), max(np.diff(zc2))])
           #if there are three traces in stream
           #get maximum of three signal periods
           if len(x) == 3:
              zc1 = []
              zc2 = []
              zc3 = []
              zc1.append(ti)
              zc2.append(ti)
              zc3.append(ti)
              i = 0
              for j in range(isignal[0][1],isignal[0][len(t[isignal]) - 1]):
                  i = i+ 1
                  if x[0].data[j-1] <= 0 and x[0].data[j] >= 0:
                     zc1.append(t[isignal][i])
                  elif x[0].data[j-1] > 0 and x[0].data[j] <= 0:
                     zc1.append(t[isignal][i])
                  if x[1].data[j-1] <= 0 and x[1].data[j] >= 0:
                     zc2.append(t[isignal][i])
                  elif x[1].data[j-1] > 0 and x[1].data[j] <= 0:
                     zc2.append(t[isignal][i])
                  if x[2].data[j-1] <= 0 and x[2].data[j] >= 0:
                     zc3.append(t[isignal][i])
                  elif x[2].data[j-1] > 0 and x[2].data[j] <= 0:
                     zc3.append(t[isignal][i])
                  if len(zc1) >= 3 and len(zc2) >= 3 and len(zc3) >= 3:
                     break
              Ts = max([max(np.diff(zc1)), max(np.diff(zc2)), max(np.diff(zc3))])

           #Ts/4 is assumed as time difference between most likely and earliest possible pick!
           self.EPick = ti - Ts/4

           if self.iplot is not None:
              plt.figure(self.iplot)
              if len(x) == 1:
                 p1, = plt.plot(t, x[0].data, 'k')
                 p2, = plt.plot(t[inoise], x[0].data[inoise])
                 p3, = plt.plot(t[isignal], x[0].data[isignal], 'r')
                 p4, = plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
                 p5, = plt.plot(zc, [0, 0, 0], '*g', markersize=14)
                 plt.legend([p1, p2, p3, p4, p5], ['Data', 'Noise Window', 'Signal Window', 'Noise Level', 'Zero Crossings'])
                 plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
                 plt.plot([ti, ti], [max(x[0].data), -max(x[0].data)], 'b', linewidth=2)
                 plt.plot([self.LPick, self.LPick], [max(x[0].data)/2, -max(x[0].data)/2], '--k')
                 plt.plot([self.EPick, self.EPick], [max(x[0].data)/2, -max(x[0].data)/2], '--k')
                 plt.xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
                 plt.yticks([])
                 plt.title('Earliest-/Latest Possible and Most Likely Pick, %s' % self.Data[0].stats.station)
              elif len(x) == 2:
                 plt.subplot(2,1,1)
                 p1, = plt.plot(t, x[0].data, 'k')
                 p2, = plt.plot(t[inoise], x[0].data[inoise])
                 p3, = plt.plot(t[isignal], x[0].data[isignal], 'r')
                 p4, = plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
                 p5, = plt.plot(zc1[0:3], [0, 0, 0], '*g', markersize=14)
                 plt.legend([p1, p2, p3, p4, p5], ['Data', 'Noise Window', 'Signal Window', 'Noise Level', 'Zero Crossings'])
                 plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
                 plt.plot([ti, ti], [max(x[0].data), -max(x[0].data)], 'b', linewidth=2)
                 plt.plot([self.LPick, self.LPick], [max(x[0].data)/2, -max(x[0].data)/2], '--k')
                 plt.plot([self.EPick, self.EPick], [max(x[0].data)/2, -max(x[0].data)/2], '--k')
                 plt.plot(zc1[0:3], [0, 0, 0], '*g')
                 plt.yticks([])
                 plt.title('Earliest-/Latest Possible and Most Likely Pick, %s' % self.Data[0].stats.station)
                 plt.subplot(2,1,2)
                 plt.plot(t, x[1].data, 'k')
                 plt.plot(t[inoise], x[1].data[inoise])
                 plt.plot(t[isignal], x[1].data[isignal], 'r')
                 plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
                 plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
                 plt.plot([ti, ti], [max(x[1].data), -max(x[1].data)], 'b', linewidth=2)
                 plt.plot([self.LPick, self.LPick], [max(x[1].data)/2, -max(x[1].data)/2], '--k')
                 plt.plot([self.EPick, self.EPick], [max(x[1].data)/2, -max(x[1].data)/2], '--k')
                 plt.plot(zc2[0:3], [0, 0, 0], '*g', markersize=14)
                 plt.xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
                 plt.yticks([])
              elif len(x) == 3:
                 plt.subplot(3,1,1)
                 p1, = plt.plot(t, x[0].data, 'k')
                 p2, = plt.plot(t[inoise], x[0].data[inoise])
                 p3, = plt.plot(t[isignal], x[0].data[isignal], 'r')
                 p4, = plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
                 p5, = plt.plot(zc1[0:3], [0, 0, 0], '*g', markersize=14)
                 plt.legend([p1, p2, p3, p4, p5], ['Data', 'Noise Window', 'Signal Window', 'Noise Level', 'Zero Crossings'])
                 plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
                 plt.plot([ti, ti], [max(x[0].data), -max(x[0].data)], 'b', linewidth=2)
                 plt.plot([self.LPick, self.LPick], [max(x[0].data)/2, -max(x[0].data)/2], '--k')
                 plt.plot([self.EPick, self.EPick], [max(x[0].data)/2, -max(x[0].data)/2], '--k')
                 plt.yticks([])
                 plt.title('Earliest-/Latest Possible and Most Likely Pick, %s' % self.Data[0].stats.station)
                 plt.subplot(3,1,2)
                 plt.plot(t, x[1].data, 'k')
                 plt.plot(t[inoise], x[1].data[inoise])
                 plt.plot(t[isignal], x[1].data[isignal], 'r')
                 plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
                 plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
                 plt.plot([ti, ti], [max(x[1].data), -max(x[1].data)], 'b', linewidth=2)
                 plt.plot([self.LPick, self.LPick], [max(x[1].data)/2, -max(x[1].data)/2], '--k')
                 plt.plot([self.EPick, self.EPick], [max(x[1].data)/2, -max(x[1].data)/2], '--k')
                 plt.plot(zc2[0:3], [0, 0, 0], '*g', markersize=14)
                 plt.yticks([])
                 plt.subplot(3,1,3)
                 plt.plot(t, x[2].data, 'k')
                 plt.plot(t[inoise], x[2].data[inoise])
                 plt.plot(t[isignal], x[2].data[isignal], 'r')
                 plt.plot([t[0], t[int(len(t)) - 1]], [nlevel, nlevel], '--k') 
                 plt.plot([t[0], t[int(len(t)) - 1]], [-nlevel, -nlevel], '--k') 
                 plt.plot([ti, ti], [max(x[2].data), -max(x[2].data)], 'b', linewidth=2)
                 plt.plot([self.LPick, self.LPick], [max(x[2].data)/2, -max(x[2].data)/2], '--k')
                 plt.plot([self.EPick, self.EPick], [max(x[2].data)/2, -max(x[2].data)/2], '--k')
                 plt.plot(zc3[0:3], [0, 0, 0], '*g', markersize=14)
                 plt.xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
                 plt.yticks([])
              plt.show()
              raw_input()
              plt.close(self.iplot)

        elif self.getpick1() == None: 
           print 'EarlLatePick: No initial onset time given! Check input!'
           return

