# -*- coding: utf-8 -*-
"""
Created Dec 2014
Implementation of the picking algorithms published and described in:

Kueperkoch, L., Meier, T., Lee, J., Friederich, W., & Egelados Working Group, 2010:
Automated determination of P-phase arrival times at regional and local distances
using higher order statistics, Geophys. J. Int., 181, 1159-1170

Kueperkoch, L., Meier, T., Bruestle, A., Lee, J., Friederich, W., & Egelados
Working Group, 2012: Automated determination of S-phase arrival times using
autoregressive prediction: application ot local and regional distances, Geophys. J. Int.,
188, 687-702.

:author: MAGS2 EP3 working group / Ludger Kueperkoch	
"""
import numpy as np
import matplotlib.pyplot as plt
from  CharFuns import *
#from pylot.core.pick.CharFuns import CharacteristicFunction
import pdb

class AutoPicking(object):
    '''
    Superclass of different, automated picking algorithms applied on a CF determined
    using AIC, HOS, or AR prediction.
    ''' 
    def __init__(self, cf, Tslope, aerr, TSNR, PickWindow, aus=None, Tsmooth=None, Pick1=None):
        '''
        :param: cf, characteristic function, on which the picking algorithm is applied
        :type: `~pylot.core.pick.CharFuns.CharacteristicFunction` object

        :param: Tslope, length of time window after pick used to determine slope 
                for quality estimation [s]
        :type: float

        :param: aerr (adjusted error), percentage of maximum of CF to determine slope for quality estimation
        :type: int

        :param: TSNR, length of time windows around pick used to determine SNR [s]
        :type: tuple (T_noise, T_gap, T_signal)

        :param: PickWindow, length of pick window [s]
        :type: float

        :param: aus ("artificial uplift of samples"), find local minimum at i if aic(i-1)*(1+aus) >= aic(i)
        :type: float

        :param: Tsmooth, length of moving smoothing window to calculate smoothed CF [s]
        :type: float

        :param: Pick1, initial (prelimenary) onset time, starting point for PragPicker
        :type: float
        '''

        assert isinstance(cf, CharacteristicFunction), "%s is not a CharacteristicFunction object" % str(cf)

        self.cf = cf.getCF()
        self.Tcf = cf.getTimeArray()
        self.dt = cf.getIncrement()
        self.setTslope(Tslope)
        self.setaerr(aerr)
        self.setTSNR(TSNR)
        self.setPickWindow(PickWindow)
        self.setaus(aus)
        self.setTsmooth(Tsmooth)
        self.setpick1(Pick1)
        self.calcPick()

    def __str__(self):
        return '''\n\t{name} object:\n
        TSlope:\t{Tslope}\n
        aerr:\t{aerr}\n
        TSNR:\t\t\t{TSNR}\n
        PickWindow:\t{PickWindow}\n
        aus:\t{aus}\n
        Tsmooth:\t{Tsmooth}\n
        Pick1:\t{Pick1}\n
        '''.format(name=type(self).__name__,
                   Tslope=self.getTslope(),
                   aerr=self.getaerr(),
                   TSNR=self.getTSNR(),
                   PickWindow=self.getPickWindow(),
                   aus=self.getaus(),
                   Tsmooth=self.getTsmooth(),
                   Pick1=self.getpick1())     

    def getTslope(self):
        return self.Tslope

    def setTslope(self, Tslope):
        self.Tslope = Tslope

    def getaerr(self):
        return self.aerr 
              
    def setaerr(self, aerr):
        self.aerr = aerr
 
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
        #taper AIC-CF to get rid off side maxima
        tap = np.hanning(len(self.cf))
        aic = tap * self.cf + max(abs(self.cf))
        #get maximum of CF as starting point
        icfmax = np.argmax(aic)
        
        #find minimum in front of maximum
        lpickwindow = int(round(self.PickWindow / self.dt))
        for i in range(icfmax - 1, max([icfmax - lpickwindow, 2]), -1): 
           if aic[i - 1] >= aic[i]:
              self.Pick = self.Tcf[i]
              break
        if self.Pick == None:
           print 'AICPicker: Could not find minimum, picking window too short?'
  
        return self.Pick

class PragPicker(AutoPicking):
    '''
    Method of pragmatic picking exploiting information given by CF.
    '''

    def calcPick(self):

        if self.getpick1() is not None:
           print 'Get onset time (pick) from HOS- or AR-CF using pragmatic picking algorithm ...'

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
           cfipick = self.cf[ipick]
           Tcfpick = self.Tcf[ipick]
           cfsmoothipick = cfsmooth[ipick]
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
           if flagpick_l > 0 and flagpick_r > 0:
              if cfpick_l <= cfpick_r:
                 self.Pick = pick_l
              else:
                 self.Pick = pick_r

        else: 
           self.Pick = None
           print 'PragPicker: No initial onset time given! Check input!'
           return

