# -*- coding: utf-8 -*-
"""
Created Dec 2014
Implementation of the picking algorithms published and described in:

Küperkoch, L., Meier, T., Lee, J., Friederich, W., & Egelados Working Group, 2010:
Automated determination of P-phase arrival times at regional and local distances
using higher order statistics, Geophys. J. Int., 181, 1159-1170

Küperkoch, L., Meier, T., Brüstle, A., Lee, J., Friederich, W., & Egelados
Working Group, 2012: Automated determination of S-phase arrival times using
autoregressive prediction: application ot local and regional distances, Geophys. J. Int.,
188, 687-702.

:author: MAGS2 EP3 working group / Ludger Küperkoch	
"""
import numpy as np
import matplotlib.pyplot as plt
import pdb

class AutoPicking(object):
    '''
    Superclass of different, automated picking algorithms applied on a CF determined
    using AIC, HOS, or AR prediction.
    ''' 
    def __init__(self, cf, Tslope, aerr, TSNR, PickWindow, peps, Tsmooth):
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

        :param: peps, find local minimum at i if aic(i-1)*(1+peps) >= aic(i)
        :type: float

        :param: Tsmooth, length of moving smoothing window to calculate smoothed CF [s]
        :type: float
        '''

        self.cf = cf.getCF()
        self.Tcf = cf.getTimeArray()
        self.dt = cf.getIncrement()
        self.setTslope(Tslope)
        self.setaerr(aerr)
        self.setTSNR(TSNR)
        self.setPickWindow(PickWindow)
        self.setpeps(peps)
        self.setTsmooth(Tsmooth)
        self.calcPick()

    def __str__(self):
        return '''\n\t{name} object:\n
        TSlope:\t{Tslope}\n
        aerr:\t{aerr}\n
        TSNR:\t\t\t{TSNR}\n
        PickWindow:\t{PickWindow}\n
        peps:\t{peps}\n
        Tsmooth:\t{Tsmooth}\n
        '''.format(name=type(self).__name__,
                   Tslope=self.getTslope(),
                   aerr=self.getaerr(),
                   TSNR=self.getTSNR(),
                   PickWindow=self.getPickWindow(),
                   peps=self.getpeps(),
                   Tsmooth=self.getTsmooth())     

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

    def getpeps(self):
        return self.peps

    def setpeps(self, peps):
        self.peps = peps

    def setTsmooth(self, Tsmooth):
        self.Tsmooth = Tsmooth

    def getTsmooth(self):
        return self.Tsmooth

    def getpick(self):
        return self.Pick

    def calcPick(self):
        self.Pick = None
 
class AICPicker(AutoPicking):
    '''
    Method to derive onset time of arriving phase based on CF
    derived from AIC.
    '''

    def calcPick(self):
 
        print 'Get onset (pick) from AIC-CF ...'

        #taper AIC-CF to get rid off side maxima
        tap = np.hanning(len(self.cf))
        aic = tap * self.cf + max(abs(self.cf))
        #get maximum of CF as starting point
        icfmax = np.argmax(aic)

        #smooth CF
        aicsmooth = np.zeros(len(aic))
        ismooth = round(self.Tsmooth / self.dt)
        if len(aic) < ismooth:
           print 'AICPicker: Tsmooth larger than AIC function!'
           self.Pick = -1
           return self.Pick
        else:
           self.Pick = -1
           for i in range(1, len(aic)):
              if i > ismooth:
                 ii1 = i - ismooth
                 aicsmooth[i] = aicsmooth[i - 1] + (aic[i] - aic[ii1]) / ismooth
              else:
                 aicsmooth[i] = np.mean(aic[0:i])

        #find common, local minimum in front of maximum
        #of smoothed and unsmoothed AIC-CF
        lpickwindow = int(round(self.PickWindow / self.dt))
        for i in range(icfmax - 1, max([icfmax - lpickwindow, 2]), -1): 
           if aic[i - 1] * (1 + self.peps) >= aic[i]:
              if aicsmooth[i - 1] * (1 + self.peps) >= aicsmooth[i]:
                 self.Pick = self.Tcf[i]
                 break

        #try again with larger peps if picking failed
        if self.Pick < 0:
           peps2 = self.peps + 0.01
           for i in range(icfmax - 1, max([icfmax - lpickwindow, 2]), -1): 
              if aic[i - 1] * (1 + peps2) >= aic[i]:
                 if aicsmooth[i - 1] * (1 + peps2) >= aicsmooth[i]:
                    self.Pick = self.Tcf[i]
                    break

class PragPicker(AutoPicking):
    '''
    Method of pragmatic picking exploiting information given by CF.
    '''

    def calcPick(self):

        print 'Get onset (pick) from HOS- or AR-CF using pragmatic picking algorithm ...'
