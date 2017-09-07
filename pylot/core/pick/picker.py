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

import warnings

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelmax
from pylot.core.pick.charfuns import CharacteristicFunction
from pylot.core.pick.utils import getnoisewin, getsignalwin


class AutoPicker(object):
    '''
    Superclass of different, automated picking algorithms applied on a CF determined
    using AIC, HOS, or AR prediction.
    '''

    warnings.simplefilter('ignore')

    def __init__(self, cf, TSNR, PickWindow, iplot=0, aus=None, Tsmooth=None, Pick1=None, fig=None):
        '''
        :param: cf, characteristic function, on which the picking algorithm is applied
        :type: `~pylot.core.pick.CharFuns.CharacteristicFunction` object

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
         EarlLatePicker
        :type: float

        '''

        assert isinstance(cf, CharacteristicFunction), "%s is not a CharacteristicFunction object" % str(cf)

        self.cf = cf.getCF()
        self.Tcf = cf.getTimeArray()
        self.Data = cf.getXCF()
        self.dt = cf.getIncrement()
        self.setTSNR(TSNR)
        self.setPickWindow(PickWindow)
        self.setiplot(iplot)
        self.setaus(aus)
        self.setTsmooth(Tsmooth)
        self.setpick1(Pick1)
        self.fig = fig
        self.calcPick()

    def __str__(self):
        return '''\n\t{name} object:\n
        TSNR:\t\t\t{TSNR}\n
        PickWindow:\t{PickWindow}\n
        aus:\t{aus}\n
        Tsmooth:\t{Tsmooth}\n
        Pick1:\t{Pick1}\n
        '''.format(name=type(self).__name__,
                   TSNR=self.getTSNR(),
                   PickWindow=self.getPickWindow(),
                   aus=self.getaus(),
                   Tsmooth=self.getTsmooth(),
                   Pick1=self.getpick1())

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

    def getSNR(self):
        return self.SNR

    def getSlope(self):
        return self.slope

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


class AICPicker(AutoPicker):
    '''
    Method to derive the onset time of an arriving phase based on CF
    derived from AIC. In order to get an impression of the quality of this inital pick,
    a quality assessment is applied based on SNR and slope determination derived from the CF,
    from which the AIC has been calculated.
    '''

    def calcPick(self):

        print('AICPicker: Get initial onset time (pick) from AIC-CF ...')

        self.Pick = None
        self.slope = None
        self.SNR = None
        plt_flag = 0
        try:
            iplot = int(self.iplot)
        except:
            if self.iplot == True or self.iplot == 'True':
               iplot = 2
            else:
               iplot = 0

        # find NaN's
        nn = np.isnan(self.cf)
        if len(nn) > 1:
            self.cf[nn] = 0
        # taper AIC-CF to get rid off side maxima
        tap = np.hanning(len(self.cf))
        aic = tap * self.cf + max(abs(self.cf))
        # smooth AIC-CF
        ismooth = int(round(self.Tsmooth / self.dt))
        aicsmooth = np.zeros(len(aic))
        if len(aic) < ismooth:
            print('AICPicker: Tsmooth larger than CF!')
            return
        else:
            for i in range(1, len(aic)):
                if i > ismooth:
                    ii1 = i - ismooth
                    aicsmooth[i] = aicsmooth[i - 1] + (aic[i] - aic[ii1]) / ismooth
                else:
                    aicsmooth[i] = np.mean(aic[1: i])
        # remove offset in AIC function
        offset = abs(min(aic) - min(aicsmooth))
        aicsmooth = aicsmooth - offset
        # get maximum of HOS/AR-CF as startimg point for searching
        # minimum in AIC function 
        icfmax = np.argmax(self.Data[0].data)

        # find minimum in AIC-CF front of maximum of HOS/AR-CF
        lpickwindow = int(round(self.PickWindow / self.dt))
        for i in range(icfmax - 1, max([icfmax - lpickwindow, 2]), -1):
            if aicsmooth[i - 1] >= aicsmooth[i]:
                self.Pick = self.Tcf[i]
                break
        # if no minimum could be found:
        # search in 1st derivative of AIC-CF
        if self.Pick is None:
            diffcf = np.diff(aicsmooth)
            # find NaN's
            nn = np.isnan(diffcf)
            if len(nn) > 1:
                diffcf[nn] = 0
            # taper CF to get rid off side maxima
            tap = np.hanning(len(diffcf))
            diffcf = tap * diffcf * max(abs(aicsmooth))
            for i in range(icfmax - 1, max([icfmax - lpickwindow, 2]), -1):
                if diffcf[i - 1] >= diffcf[i]:
                    self.Pick = self.Tcf[i]
                    break

        # quality assessment using SNR and slope from CF
        if self.Pick is not None:
            # get noise window
            inoise = getnoisewin(self.Tcf, self.Pick, self.TSNR[0], self.TSNR[1])
            # check, if these are counts or m/s, important for slope estimation!
            # this is quick and dirty, better solution?
            if max(self.Data[0].data < 1e-3) and max(self.Data[0].data >= 1e-6):
                self.Data[0].data = self.Data[0].data * 1000000.
            elif max(self.Data[0].data < 1e-6):
                self.Data[0].data = self.Data[0].data * 1e13
            # get signal window
            isignal = getsignalwin(self.Tcf, self.Pick, self.TSNR[2])
            if len(isignal) == 0:
                return
            ii = min([isignal[len(isignal) - 1], len(self.Tcf)])
            isignal = isignal[0:ii]
            try:
                self.Data[0].data[isignal]
            except IndexError as e:
                msg = "Time series out of bounds! {}".format(e)
                print(msg)
                return
            # calculate SNR from CF
            self.SNR = max(abs(self.Data[0].data[isignal] - np.mean(self.Data[0].data[isignal]))) / \
                       max(abs(self.Data[0].data[inoise] - np.mean(self.Data[0].data[inoise])))
            # calculate slope from CF after initial pick
            # get slope window
            tslope = self.TSNR[3]  # slope determination window
            islope = np.where((self.Tcf <= min([self.Pick + tslope, self.Tcf[-1]])) \
                              & (self.Tcf >= self.Pick)) # TODO: put this in a seperate function like getsignalwin
            # find maximum within slope determination window
            # 'cause slope should be calculated up to first local minimum only!
            try:
                dataslope = self.Data[0].data[islope[0][0:-1]]
            except IndexError:
                print("Slope Calculation: empty array islope, check signal window")
                return
            if len(dataslope) < 1:
                print('No data in slope window found!')
                return
            imaxs, = argrelmax(dataslope)
            if imaxs.size:
                imax = imaxs[0]
            else:
                imax = np.argmax(dataslope)
            iislope = islope[0][0:imax + 1]
            if len(iislope) < 2:
                # calculate slope from initial onset to maximum of AIC function
                print("AICPicker: Not enough data samples left for slope calculation!")
                print("Calculating slope from initial onset to maximum of AIC function ...")
                imax = np.argmax(aicsmooth[islope[0][0:-1]])
                if imax == 0:
                    print("AICPicker: Maximum for slope determination right at the beginning of the window!")
                    print("Choose longer slope determination window!")
                    if self.iplot > 1:
                        if self.fig == None or self.fig == 'None':
                            fig = plt.figure()  # self.iplot) ### WHY? MP MP
                            plt_flag = 1
                        else:
                            fig = self.fig
                        ax = fig.add_subplot(111)
                        x = self.Data[0].data
                        ax.plot(self.Tcf, x / max(x), 'k', label='(HOS-/AR-) Data')
                        ax.plot(self.Tcf, aicsmooth / max(aicsmooth), 'r', label='Smoothed AIC-CF')
                        ax.legend(loc=1)
                        ax.set_xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
                        ax.set_yticks([])
                        ax.set_title(self.Data[0].stats.station)
                        if plt_flag == 1:
                            fig.show()
                            try: input()
                            except SyntaxError: pass
                            plt.close(fig)
                    return
                iislope = islope[0][0:imax+1]
            dataslope = self.Data[0].data[iislope]
            # calculate slope as polynomal fit of order 1
            xslope = np.arange(0, len(dataslope), 1)
            P = np.polyfit(xslope, dataslope, 1)
            datafit = np.polyval(P, xslope)
            if datafit[0] >= datafit[-1]:
                print('AICPicker: Negative slope, bad onset skipped!')
                return
            self.slope = 1 / (len(dataslope) * self.Data[0].stats.delta) * (datafit[-1] - datafit[0])

        else:
            self.SNR = None
            self.slope = None

        if iplot > 1:
            if self.fig == None or self.fig == 'None':
                fig = plt.figure()  # self.iplot)
                plt_flag = 1
            else:
                fig = self.fig
            ax1 = fig.add_subplot(211)
            x = self.Data[0].data
            if len(self.Tcf) > len(self.Data[0].data): # why? LK
                self.Tcf = self.Tcf[0:len(self.Tcf)-1]
            ax1.plot(self.Tcf, x / max(x), 'k', label='(HOS-/AR-) Data')
            ax1.plot(self.Tcf, aicsmooth / max(aicsmooth), 'r', label='Smoothed AIC-CF')
            if self.Pick is not None:
                ax1.plot([self.Pick, self.Pick], [-0.1, 0.5], 'b', linewidth=2, label='AIC-Pick')
            ax1.set_xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
            ax1.set_yticks([])
            ax1.legend(loc=1)

            if self.Pick is not None:
                ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)
                ax2.plot(self.Tcf, x, 'k', label='Data')
                ax1.axvspan(self.Tcf[inoise[0]], self.Tcf[inoise[-1]], color='y', alpha=0.2, lw=0, label='Noise Window')
                ax1.axvspan(self.Tcf[isignal[0]], self.Tcf[isignal[-1]], color='b', alpha=0.2, lw=0,
                            label='Signal Window')
                ax1.axvspan(self.Tcf[iislope[0]], self.Tcf[iislope[-1]], color='g', alpha=0.2, lw=0,
                            label='Slope Window')

                ax2.axvspan(self.Tcf[inoise[0]], self.Tcf[inoise[-1]], color='y', alpha=0.2, lw=0, label='Noise Window')
                ax2.axvspan(self.Tcf[isignal[0]], self.Tcf[isignal[-1]], color='b', alpha=0.2, lw=0,
                            label='Signal Window')
                ax2.axvspan(self.Tcf[iislope[0]], self.Tcf[iislope[-1]], color='g', alpha=0.2, lw=0,
                            label='Slope Window')
                ax2.plot(self.Tcf[iislope], datafit, 'g', linewidth=2, label='Slope')

                ax1.set_title('Station %s, SNR=%7.2f, Slope= %12.2f counts/s' % (self.Data[0].stats.station,
                                                                                 self.SNR, self.slope))
                ax2.set_xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
                ax2.set_ylabel('Counts')
                ax2.set_yticks([])
                ax2.legend(loc=1)
                if plt_flag == 1:
                    fig.show()
                    try: input()
                    except SyntaxError: pass
                    plt.close(fig)
            else:
                ax1.set_title(self.Data[0].stats.station)
                if plt_flag == 1:
                    fig.show()
                    try: input()
                    except SyntaxError: pass
                    plt.close(fig)

        if self.Pick == None:
            print('AICPicker: Could not find minimum, picking window too short?')

        return


class PragPicker(AutoPicker):
    '''
    Method of pragmatic picking exploiting information given by CF.
    '''

    def calcPick(self):

        try:
            iplot = int(self.getiplot())
        except:
            if self.getiplot() == True or self.getiplot() == 'True':
               iplot = 2
            else:
               iplot = 0

        if self.getpick1() is not None:
            print('PragPicker: Get most likely pick from HOS- or AR-CF using pragmatic picking algorithm ...')

            self.Pick = None
            self.SNR = None
            self.slope = None
            pickflag = 0
            plt_flag = 0
            # smooth CF
            ismooth = int(round(self.Tsmooth / self.dt))
            cfsmooth = np.zeros(len(self.cf))
            if len(self.cf) < ismooth:
                print('PragPicker: Tsmooth larger than CF!')
                return
            else:
                for i in range(1, len(self.cf)):
                    if i > ismooth:
                        ii1 = i - ismooth
                        cfsmooth[i] = cfsmooth[i - 1] + (self.cf[i] - self.cf[ii1]) / ismooth
                    else:
                        cfsmooth[i] = np.mean(self.cf[1: i])

            # select picking window
            # which is centered around tpick1
            ipick = np.where((self.Tcf >= self.getpick1() - self.PickWindow / 2) \
                             & (self.Tcf <= self.getpick1() + self.PickWindow / 2))
            cfipick = self.cf[ipick] - np.mean(self.cf[ipick])
            Tcfpick = self.Tcf[ipick]
            cfsmoothipick = cfsmooth[ipick] - np.mean(self.cf[ipick])
            ipick1 = np.argmin(abs(self.Tcf - self.getpick1()))
            cfpick1 = 2 * self.cf[ipick1]

            # check trend of CF, i.e. differences of CF and adjust aus ("artificial uplift 
            # of picks") regarding this trend
            # prominent trend: decrease aus
            # flat: use given aus
            cfdiff = np.diff(cfipick)
            if len(cfdiff)<20:
                print('PragPicker: Very few samples for CF. Check LTA window dimensions!')
            i0diff = np.where(cfdiff > 0)
            cfdiff = cfdiff[i0diff]
            if len(cfdiff)<1:
                print('PragPicker: Negative slope for CF. Check LTA window dimensions! STOP')
                self.Pick = None
                return
            minaus = min(cfdiff * (1 + self.aus))
            aus1 = max([minaus, self.aus])

            # at first we look to the right until the end of the pick window is reached
            flagpick_r = 0
            flagpick_l = 0
            cfpick_r = 0
            cfpick_l = 0
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

            # now we look to the left
            if len(self.cf) > ipick1 +1:
                for i in range(ipick1, max([ipick1 - lpickwindow + 1, 2]), -1):
                    if self.cf[i + 1] > self.cf[i] and self.cf[i - 1] >= self.cf[i]:
                        if cfsmooth[i - 1] * (1 + aus1) >= cfsmooth[i]:
                            if cfpick1 >= self.cf[i]:
                                pick_l = self.Tcf[i]
                                self.Pick = pick_l
                                flagpick_r = 1
                                cfpick_l = self.cf[i]
                                break
            else:
                msg ='PragPicker: Initial onset too close to start of CF! \
                      Stop finalizing pick to the left.'
                print(msg)

            # now decide which pick: left or right?
            if flagpick_l > 0 and flagpick_r > 0 and cfpick_l <= 3 * cfpick_r:
                self.Pick = pick_l
                pickflag = 1
            elif flagpick_l > 0 and flagpick_r > 0 and cfpick_l >= cfpick_r:
                self.Pick = pick_r
                pickflag = 1
            elif flagpick_l == 0 and flagpick_r > 0 and cfpick_l >= cfpick_r:
                self.Pick = pick_l
                pickflag = 1
            else:
                print("PragPicker: Could not find reliable onset!")
                self.Pick = None
                pickflag = 0

            if iplot > 1:
                if self.fig == None or self.fig == 'None':
                    fig = plt.figure()  # self.getiplot())
                    plt_flag = 1
                else:
                    fig = self.fig
                ax = fig.add_subplot(111)
                ax.plot(Tcfpick, cfipick, 'k', label='CF')
                ax.plot(Tcfpick, cfsmoothipick, 'r', label='Smoothed CF')
                if pickflag > 0:
                    ax.plot([self.Pick, self.Pick], [min(cfipick), max(cfipick)], 'b', linewidth=2, label='Pick')
                ax.set_xlabel('Time [s] since %s' % self.Data[0].stats.starttime)
                ax.set_yticks([])
                ax.set_title(self.Data[0].stats.station)
                ax.legend(loc=1)
                if plt_flag == 1:
                    fig.show()
                    try: input()
                    except SyntaxError: pass
                    plt.close(fig)
                return

        else:
            print("PragPicker: No initial onset time given! Check input!")
            self.Pick = None
            return
