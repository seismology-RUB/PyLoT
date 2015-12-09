#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
from obspy.core import read
from obspy import Stream
from obspy import Trace
from pylot.core.pick.CharFuns import HOScf
from pylot.core.pick.CharFuns import AICcf
from pylot.core.pick.utils import getSNR
from pylot.core.pick.utils import earllatepicker
import matplotlib.pyplot as plt
plt.interactive('True')

class SeismicShot(object):
    '''
    SuperClass for a seismic shot object.
    '''
    def __init__(self, obsfile):
        '''
        Initialize seismic shot object giving an inputfile.

        :param: obsfile, ((!SEG2/SEGY!)) file readable by obspy
        :type: string
        '''
        self.traces = read(obsfile)
        self.recCoordlist = None
        self.srcCoordlist = None
        self.traceIDs = None
        self.picks = {}
        self.pwindow= {}
        self.manualpicks= {}
        self.snr = {}
        self.snrthreshold = {}
        self.timeArray = {}
        self.traces4plot = {}
        self.paras = {}
        self.paras['shotname'] = obsfile

    def removeEmptyTraces(self):
        traceIDs = []
        coordlist = self.getRecCoordlist()
        removed = []
        for i in range(0, len(coordlist)):
            traceIDs.append(int(coordlist[i].split()[0]))

        for trace in self.traces:
            try:
                traceIDs.index(int(trace.stats.channel))
            except:
                self.traces.remove(trace)
                removed.append(int(trace.stats.channel))

        if len(removed) > 0:
            return removed

    def removeTrace(self, traceID):
        for trace in self.traces:
            if traceID == trace.stats.channel:
                self.traces.remove(trace)

        # for traceID in TraceIDs:
        #     traces = [trace for trace in self.traces if int(trace.stats.channel) == traceID]
        #     if len(traces) is not 1:
        #         self.traces.remove(trace)

    def updateTraceList(self):
        '''
        Looks for empty traces, returns a list of deleted traceIDs.
        '''
        traceIDs = []
        for traceID in self.getTraceIDlist():
            if traceID not in self.getStreamTraceIDs():
                self.traceIDs.remove(traceID)
                traceIDs.append(traceID)
        return traceIDs

    def setParameters(self, name, value):
        self.paras[name] = value

    def setCut(self, cut):
        self.setParameters('cut', cut)

    def setTmovwind(self, tmovwind):
        self.setParameters('tmovwind', tmovwind)

    def setOrder(self, order):
       self.setParameters('order', order)

    def setTsignal(self, tsignal):
       self.setParameters('tsignal', tsignal)

    def setTgap(self, tgap):
       self.setParameters('tgap', tgap)

    def setShotnumber(self, shotname):
       self.setParameters('shotname', shotname)

    def setRecfile(self, recfile):
       self.setParameters('recfile', recfile)

    def setSourcefile(self, sourcefile):
       self.setParameters('sourcefile', sourcefile)

    def getParas(self):
        return self.paras

    def getShotname(self):
        return self.paras['shotname']

    def getCut(self):
        return self.paras['cut']

    def getTmovwind(self):
        return self.paras['tmovwind']

    def getOrder(self):
        return self.paras['order']

    def getTsignal(self):
        return self.paras['tsignal']

    def getTgap(self):
        return self.paras['tgap']

    def getShotnumber(self):
        return self.paras['shotnumber']

    def getRecfile(self):
        return self.paras['recfile']

    def getSourcefile(self):
        return self.paras['sourcefile']

    def getManualPick(self, traceID):
        if not self.getManualPickFlag(traceID) == 0:
            return self.manualpicks[traceID]['mpp']

    def getManualEarliest(self, traceID):
        return self.manualpicks[traceID]['epp']

    def getManualLatest(self, traceID):
        return self.manualpicks[traceID]['lpp']

    def getPick(self, traceID, returnRemoved = False):
        if not self.getPickFlag(traceID) == 0:
            return self.picks[traceID]['mpp']
        if returnRemoved == True:
            #print('getPick: Returned removed pick for shot %d, traceID %d' %(self.getShotnumber(), traceID))
            return self.picks[traceID]['mpp']

    def getPickIncludeRemoved(self, traceID):
        return self.getPick(traceID, returnRemoved = True)

    def getEarliest(self, traceID):
        return self.picks[traceID]['epp']

    def getLatest(self, traceID):
        return self.picks[traceID]['lpp']

    def getSymmetricPickError(self, traceID):
        pickerror = self.picks[traceID]['spe']
        if np.isnan(pickerror) == True:
            print "SPE is NaN for shot %s, traceID %s"%(self.getShotnumber(), traceID)
        return pickerror

    def getPickError(self, traceID):
        pickerror = abs(self.getEarliest(traceID) - self.getLatest(traceID)) / 2
        if np.isnan(pickerror) == True:
            print("PE is NaN for shot %s, traceID %s"%(self.getShotnumber(), traceID))
        return pickerror

    def getStreamTraceIDs(self):
        traceIDs = []
        for trace in self.traces:
            traceIDs.append(int(trace.stats.channel))
        return traceIDs

    def getTraceIDlist(self):
        '''
        Returns a list containing the traceIDs read from the receiver inputfile.
        '''
        traceIDs = []
        if self.traceIDs == None:
            recCoordlist = self.getRecCoordlist()
            for i in range(0, len(recCoordlist)):
                traceIDs.append(int(recCoordlist[i].split()[0]))
            self.traceIDs = traceIDs

        return self.traceIDs

    def getPickwindow(self, traceID):
        try:
            self.pwindow[traceID]
        except KeyError as e:
            print('no pickwindow for trace %s, set to %s' % (traceID, self.getCut()))
            self.setPickwindow(traceID, self.getCut())
        return self.pwindow[traceID]

    def getSNR(self, traceID):
        return self.snr[traceID]

    def getSNRthreshold(self, traceID):
        return self.snrthreshold[traceID]

    def getRecCoordlist(self):
        if self.recCoordlist is None:
            coordlist = open(self.getRecfile(),'r').readlines()
            #print 'Reading receiver coordinates from %s' %(self.getRecfile())
            self.recCoordlist = coordlist
        return self.recCoordlist

    def getSrcCoordlist(self):
        if self.srcCoordlist is None:
            coordlist = open(self.getSourcefile(),'r').readlines()
            #print 'Reading shot coordinates from %s' %(self.getSourcefile())
            self.srcCoordlist = coordlist
        return self.srcCoordlist

    def getTimeArray(self, traceID):
        return self.timeArray[traceID]

    def getHOScf(self, traceID):
        '''
        Returns the higher order statistics characteristic function for a trace using pylot.

        :param: traceID
        :type: int

        :param: cut, cut out a part of the trace (t_start, t_end) [s]
        :type: tuple

        :param: t2, size of the moving window [s]
        :type: float

        :param: order, order of the characteristic function
        :type: int
        '''
        return HOScf(self.getSingleStream(traceID), self.getCut(),
                     self.getTmovwind(), self.getOrder(), stealthMode = True)

    def getAICcf(self, traceID):
        '''
        Returns the Akaike criterion for a trace using pylot and the higher order statistics characteristic function.

        :param: traceID
        :type: int

        :param: cut, cut out a part of the trace (t_start, t_end) [s]
        :type: tuple

        :param: t2, size of the moving window [s]
        :type: float

        :param: order, order of the characteristic function
        :type: int
        '''

        st_cf = Stream()
        tr_cf = Trace()
        tr_cf.data = self.getHOScf(traceID).getCF()
        st_cf += tr_cf
        return AICcf(st_cf, self.getCut(), self.getTmovwind(), stealthMode = True)

    def getSingleStream(self, traceID):  ########## SEG2 / SEGY ? ##########
        '''
        Returns a Stream with only one trace (instead of just one trace) because it is needed for pylot.

        :param: traceID
        :type: int
        '''
        #traces = [trace for trace in self.traces if int(trace.stats.seg2['CHANNEL_NUMBER']) == traceID]
        traces = [trace for trace in self.traces if int(trace.stats.channel) == traceID]
        if len(traces) == 1:
            return Stream(traces)
        self.setPick(traceID, None)
        print 'Warning: ambigious or empty traceID: %s' % traceID

        #raise ValueError('ambigious or empty traceID: %s' % traceID)

    def pickTraces(self, traceID, windowsize, folm = 0.6, HosAic = 'hos'): ########## input variables ##########
        # LOCALMAX NOT IMPLEMENTED!
        '''
        Intitiate picking for a trace.

        :param: traceID
        :type: int

        :param: cutwindow (equals HOScf 'cut' variable)
        :type: tuple

        :param: t2 (equals HOScf t2 variable)
        :type: float

        :param: order (equals HOScf 'order' variable)
        :type: int

        :param: windowsize, window around the returned HOS picktime, to search for the AIC minumum
        :type: 'tuple'

        :param: folm, fraction of local maximumm (default = 0.6)
        :type: 'real'

        :param: HosAic, get hos or aic pick (can be 'hos'(default) or 'aic')
        :type: 'string'
        '''
        hoscf = self.getHOScf(traceID) ### determination of both, HOS and AIC (need to change threshold-picker) ###
        aiccf = self.getAICcf(traceID)

        self.timeArray[traceID] = hoscf.getTimeArray()
        aiccftime, hoscftime = self.threshold(hoscf, aiccf, windowsize, self.getPickwindow(traceID), folm)
        setHosAic = {'hos': hoscftime,
                     'aic': aiccftime}

        self.setPick(traceID, setHosAic[HosAic])

    def setEarllatepick(self, traceID, nfac = 1.5):
        tgap = self.getTgap()
        tsignal = self.getTsignal()
        tnoise = self.getPickIncludeRemoved(traceID) - tgap

        (self.picks[traceID]['epp'],
         self.picks[traceID]['lpp'],
         self.picks[traceID]['spe']) = earllatepicker(self.getSingleStream(traceID),
                                                     nfac, (tnoise, tgap, tsignal),
                                                     self.getPickIncludeRemoved(traceID),
                                                     stealthMode = True)

        if self.picks[traceID]['epp'] < 0:
            self.picks[traceID]['epp']
            #print('setEarllatepick: Set epp to 0 because it was < 0')

        # TEST OF 1/2 PICKERROR
        # self.picks[traceID]['spe'] *= 0.5
        # TEST OF 1/2 PICKERROR

    def threshold(self, hoscf, aiccf, windowsize, pickwindow, folm = 0.6):
        '''
        Threshold picker, using the local maximum in a pickwindow to find the time at
        which a fraction of the local maximum is reached for the first time.

        :param: hoscf, Higher Order Statistics Characteristic Function
        :type: 'Characteristic Function'

        :param: aiccf, Characteristic Function after Akaike
        :type: 'Characteristic Function'

        :param: windowsize, window around the returned HOS picktime, to search for the AIC minumum
        :type: 'tuple'

        :param: pickwindow [seconds]
        :type: 'tuple'

        :param: cutwindow [seconds], cut a part of the trace as in Characteristic Function
        :type: 'tuple'

        :param: folm, fraction of local maximum (default = 0.6)
        :type: 'real'
        '''
        hoscflist = list(hoscf.getCF())
        leftb = int(pickwindow[0] / self.getCut()[1] * len(hoscflist))
        rightb = int(pickwindow[1] / self.getCut()[1] * len(hoscflist))

        threshold = folm * max(hoscflist[leftb : rightb]) # combination of local maximum and threshold

        m = leftb

        while hoscflist[m] < threshold:
            m += 1

        hoscftime = list(hoscf.getTimeArray())[m]

        lb = max(0, m - windowsize[0]) # if window exceeds t = 0
        aiccfcut = list(aiccf.getCF())[lb : m + windowsize[1]]
        n = aiccfcut.index(min(aiccfcut))

        m = lb + n

        aiccftime = list(hoscf.getTimeArray())[m]

        return aiccftime, hoscftime

    def getDistance(self, traceID):
        '''
        Returns the distance of the receiver with the ID == traceID to the source location (shot location).
        Uses getSrcLoc and getRecLoc.

        :param: traceID
        :type: int
        '''
        shotX, shotY, shotZ = self.getSrcLoc()
        recX, recY, recZ = self.getRecLoc(traceID)
        dist = np.sqrt((shotX-recX)**2 + (shotY-recY)**2 + (shotZ-recZ)**2)

        if np.isnan(dist) == True:
            raise ValueError("Distance is NaN for traceID %s" %traceID)

        return dist
        #return abs(float(self.getSrcLoc(traceID))-float(self.getRecLoc(traceID)))

    def getRecLoc(self, traceID):  ########## input FILENAME ##########
        '''
        Returns the location (x, y, z) of the receiver with the ID == traceID.
        RECEIVEIVER FILE MUST BE SET FIRST, TO BE IMPROVED.

        :param: traceID
        :type: int
        '''
        if traceID == 0:        # artificial traceID 0 with pick at t = 0
            return self.getSrcLoc()

        coordlist = self.getRecCoordlist()
        for i in range(0, len(coordlist)):
            if int(coordlist[i].split()[0]) == traceID:
                x = coordlist[i].split()[1]
                y = coordlist[i].split()[2]
                z = coordlist[i].split()[3]
                return float(x), float(y), float(z)

        #print "WARNING: traceID %s not found" % traceID
        raise ValueError("traceID %s not found" % traceID)
        #return float(self.getSingleStream(traceID)[0].stats.seg2['RECEIVER_LOCATION'])

    def getSrcLoc(self):  ########## input FILENAME ##########
        '''
        Returns the location (x, y, z) of the shot.
        SOURCE FILE MUST BE SET FIRST, TO BE IMPROVED.
        '''
        coordlist = self.getSrcCoordlist()
        for i in range(0, len(coordlist)):
            if int(coordlist[i].split()[0]) == self.paras['shotnumber']:
                x = coordlist[i].split()[1]
                y = coordlist[i].split()[2]
                z = coordlist[i].split()[3]
                return float(x), float(y), float(z)
        #return float(self.getSingleStream(traceID)[0].stats.seg2['SOURCE_LOCATION'])

    def getTraceIDs4Dist(self, distance = 0, distancebin = (0, 0)): ########## nur fuer 2D benutzt, 'distance bins' ##########
        '''
        Returns the traceID(s) for a certain distance between source and receiver.
        Used for 2D Tomography. TO BE IMPROVED.

        :param: distance
        :type: real

        :param: distancebin
        :type: tuple
        '''

        traceID_list = []
        for trace in self.traces:
                #traceID = int(trace.stats.seg2['CHANNEL_NUMBER'])
                traceID = int(trace.stats.channel)
                if distance != 0:
                    if self.getDistance(traceID) == distance:
                        traceID_list.append(traceID)
                if distancebin[0] >= 0 and distancebin[1] > 0:
                    if distancebin[0] < self.getDistance(traceID) < distancebin[1]:
                        traceID_list.append(traceID)

        if len(traceID_list) > 0:
            return traceID_list

    # def setManualPicks(self, traceID, picklist): ########## picklist momentan nicht allgemein, nur testweise benutzt ##########
    #     '''
    #     Sets the manual picks for a receiver with the ID == traceID for comparison.

    #     :param: traceID
    #     :type: int

    #     :param: picklist, list containing the manual picks (mostlikely, earliest, latest).
    #     :type: list
    #     '''
    #     picks = picklist[traceID - 1].split()
    #     mostlikely = float(picks[1])
    #     earliest = float(picks[2])
    #     latest = float(picks[3])

    #     if not self.manualpicks.has_key(traceID):
    #         self.manualpicks[traceID] = (mostlikely, earliest, latest)
        #else:
        #    raise KeyError('MANUAL pick to be set more than once for traceID %s' % traceID)

    def setManualPicksFromFile(self, directory = 'picks'):
        '''
        Read manual picks from *.pck file.
        The * must be identical with the shotnumber.
        '''
        if directory[-1] == '/':
            filename = directory + str(self.getShotnumber()) + '.pck'
        else:
            filename = directory + '/' + str(self.getShotnumber()) + '.pck'
        infile = open(filename, 'r')
        mpicks = infile.readlines()
        for line in mpicks:
            if line.split()[0] == []:
                continue
            traceID, mpp, epp, lpp = line.split()
            traceID = int(traceID)
            if traceID in self.picks.keys():
                self.manualpicks[traceID] = {'mpp': float(mpp),
                                             'epp': float(epp),
                                             'lpp': float(lpp)}
            if float(mpp) <= 0:
                self.setManualPickFlag(traceID, 0)
            else:
                self.setManualPickFlag(traceID, 1)


    def setPick(self, traceID, pick): ########## siehe Kommentar ##########
        if not traceID in self.picks.keys():
            self.picks[traceID] = {}
        self.picks[traceID]['mpp'] = pick
        self.picks[traceID]['flag'] = 1
        # ++++++++++++++ Block raus genommen, da Error beim 2ten Mal picken! (Ueberschreiben von erstem Pick!)
        # if not self.picks.has_key(traceID):
        #     self.getPick(traceID) = picks
        # else:
        #     raise KeyError('pick to be set more than once for traceID %s' % traceID)

        #    def readParameter(self, parfile):
        #        parlist = open(parfile,'r').readlines()

    def removePick(self, traceID):
        self.setPickFlag(traceID, 0)

    def setPickFlag(self, traceID, flag):
        'Set flag = 0 if pick is invalid, else flag = 1'
        self.picks[traceID]['flag'] = flag

    def getPickFlag(self, traceID):
        return self.picks[traceID]['flag']

    def setManualPickFlag(self, traceID, flag):
        'Set flag = 0 if pick is invalid, else flag = 1'
        self.manualpicks[traceID]['flag'] = flag

    def getManualPickFlag(self, traceID):
        return self.manualpicks[traceID]['flag']

    def setPickwindow(self, traceID, pickwindow):
        self.pwindow[traceID] = pickwindow

    def setSNR(self, traceID):  ########## FORCED HOS PICK ##########
        '''
        Gets the SNR using pylot and then sets the SNR for the traceID.

        :param: traceID
        :type: int

        :param: (tnoise, tgap, tsignal), as used in pylot SNR
        '''

        from pylot.core.pick.utils import getSNR

        tgap = self.getTgap()
        tsignal = self.getTsignal()
        tnoise = self.getPick(traceID) - tgap

        self.snr[traceID] = getSNR(self.getSingleStream(traceID), (tnoise,tgap,tsignal), self.getPick(traceID))

    def setSNRthreshold(self, traceID, snrthreshold):
        self.snrthreshold[traceID] = snrthreshold

    def getDistArray4ttcPlot(self):  ########## nur fuer 2D benoetigt ##########
        '''
        Function to create a distance array for the plots. 2D only! X DIRECTION!!
        '''
        distancearray = []

        for traceID in self.picks.keys():
            if self.getRecLoc(traceID)[0] > self.getSrcLoc()[0]:
                distancearray.append(self.getDistance(traceID))
            elif self.getRecLoc(traceID)[0] <= self.getSrcLoc()[0]:
                distancearray.append((-1)*self.getDistance(traceID))

        return distancearray


    def plot2dttc(self, ax = None):  ########## 2D ##########
        '''
        Function to plot the traveltime curve for automated picks of a shot. 2d only! ATM: X DIRECTION!!
        '''
        import matplotlib.pyplot as plt
        plt.interactive('True')
        picks = []

        for traceID in self.picks.keys():
            picks.append(self.getPick(traceID))

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        # shotnumbers = [shotnumbers for (shotnumbers, shotnames) in sorted(zip(shotnumbers, shotnames))]
        plotarray = sorted(zip(self.getDistArray4ttcPlot(), picks))
        x = []; y = []
        for point in plotarray:
            x.append(point[0])
            y.append(point[1])
        ax.plot(x, y,'r', label = "Automatic Picks")
        plt.title('Shot: %s' %self.getShotnumber())

    def plotmanual2dttc(self, ax = None):  ########## 2D ##########
        '''
        Function to plot the traveltime curve for manual picks of a shot. 2D only!
        '''
        import matplotlib.pyplot as plt
        plt.interactive('True')
        manualpicktimesarray = []

        for traceID in self.picks.keys():
            if not traceID in self.manualpicks.keys() or self.getManualPickFlag(traceID) == 0:
                manualpicktimesarray.append(None)
            else:
                manualpicktimesarray.append(self.getManualPick(traceID))

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        plotarray = sorted(zip(self.getDistArray4ttcPlot(), manualpicktimesarray))
        x = []; y = []
        for point in plotarray:
            x.append(point[0])
            y.append(point[1])
        ax.plot(x, y, 'b', label = "Manual Picks")

    # def plotpickwindow(self):   ########## 2D ##########
    #     '''
    #     Plots the pickwindow of a shot for the 2nd iteration step of picking. 2D only!
    #     '''
    #     import matplotlib.pyplot as plt
    #     plt.interactive('True')
    #     pickwindowarray_lowerb = []
    #     pickwindowarray_upperb = []

    #     i = 1
    #     for traceID in self.pwindow.keys():
    #         pickwindowarray_lowerb.append(self.pwindow[traceID][0])
    #         pickwindowarray_upperb.append(self.pwindow[traceID][1])
    #         i += 1

    #     plt.plot(self.getDistArray4ttcPlot(), pickwindowarray_lowerb, ':k')
    #     plt.plot(self.getDistArray4ttcPlot(), pickwindowarray_upperb, ':k')

    def plotTrace(self, traceID, plotSNR = True, lw = 1):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax = self._drawStream(traceID, ax = ax)

        tgap = self.getTgap()
        tsignal = self.getTsignal()
        pick = self.getPick(traceID)
        tnoise = pick - tgap
        snr, snrdb, noiselevel = self.getSNR(traceID)

        ax.plot([0, tnoise], [noiselevel, noiselevel], 'm', linewidth = lw, label = 'noise level')
        ax.plot([tnoise, pick], [noiselevel, noiselevel], 'g:', linewidth = lw, label = 'gap')
        ax.plot([tnoise + tgap, pick + tsignal], [noiselevel * snr, noiselevel * snr], 'b', linewidth = lw, label = 'signal level')
        ax.legend()
        ax.text(0.05, 0.95, 'SNR: %s' %snr, transform = ax.transAxes)

    def plot_traces(self, traceID, folm = 0.6): ########## 2D, muss noch mehr verbessert werden ##########
        from matplotlib.widgets import Button

        def onclick(event):
            self.setPick(traceID, event.xdata)
            if self.getSNR(traceID)[0] > 1:
                self.setEarllatepick(traceID)
            self._drawStream(traceID, refresh = True)
            self._drawCFs(traceID, folm, refresh = True)
            fig.canvas.mpl_disconnect(self.traces4plot[traceID]['cid'])
            plt.draw()

        def rmPick(event = None):
            self.removePick(traceID)
            self._drawStream(traceID, refresh = True)
            self._drawCFs(traceID, folm, refresh = True)
            plt.draw()

        def connectButton(event = None):
            cid = fig.canvas.mpl_connect('button_press_event', onclick)
            self.traces4plot[traceID]['cid'] = cid

        def cleanup(event):
            self.traces4plot[traceID] = {}

        fig = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2, sharex = ax1)
        axb1 = fig.add_axes([0.15, 0.91, 0.05, 0.03])
        axb2 = fig.add_axes([0.22, 0.91, 0.05, 0.03])
        button1 = Button(axb1, 'repick', color = 'red', hovercolor = 'grey')
        button1.on_clicked(connectButton)
        button2 = Button(axb2, 'delete', color = 'green', hovercolor = 'grey')
        button2.on_clicked(rmPick)
        fig.canvas.mpl_connect('close_event', cleanup)

        self.traces4plot[traceID] = dict(fig=fig, ax1=ax1, ax2=ax2, axb1=axb1, axb2=axb2, button1=button1,
                                         button2=button2, cid=None)

        self._drawStream(traceID)
        self._drawCFs(traceID, folm)

    def _drawStream(self, traceID, refresh = False, ax = None):
        from pylot.core.util.utils import getGlobalTimes
        from pylot.core.util.utils import prepTimeAxis

        stream = self.getSingleStream(traceID)
        stime = getGlobalTimes(stream)[0]
        timeaxis = prepTimeAxis(stime, stream[0])
        timeaxis -= stime

        if ax is None:
            ax = self.traces4plot[traceID]['ax1']

        if refresh == True:
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.clear()
        if refresh == True:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        ax.set_title('Shot: %s, traceID: %s, pick: %s'
                     %(self.getShotnumber(), traceID, self.getPick(traceID)))
        ax.plot(timeaxis, stream[0].data, 'k', label = 'trace')
        ax.plot([self.getPick(traceID), self.getPick(traceID)],
                [ax.get_ylim()[0],
                 ax.get_ylim()[1]],
                'r', label = 'most likely')
        if self.getEarliest(traceID) is not None:
            ax.plot([self.getEarliest(traceID), self.getEarliest(traceID)],
                    [ax.get_ylim()[0],
                     ax.get_ylim()[1]],
                    'g:', label = 'earliest')
        if self.getLatest(traceID) is not None:
            ax.plot([self.getLatest(traceID), self.getLatest(traceID)],
                    [ax.get_ylim()[0],
                     ax.get_ylim()[1]],
                    'b:', label = 'latest')

        ax.legend()
        return ax

    def _drawCFs(self, traceID, folm, refresh = False):
        hoscf = self.getHOScf(traceID)
        aiccf = self.getAICcf(traceID)
        ax = self.traces4plot[traceID]['ax2']

        if refresh == True:
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.clear()
        if refresh == True:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        ax.plot(hoscf.getTimeArray(), hoscf.getCF(), 'b', label = 'HOS')
        ax.plot(hoscf.getTimeArray(), aiccf.getCF(), 'g', label = 'AIC')
        ax.plot([self.getPick(traceID), self.getPick(traceID)],
                [ax.get_ylim()[0],
                 ax.get_ylim()[1]],
                'r', label = 'most likely')
        if self.getEarliest(traceID) is not None:
            ax.plot([self.getEarliest(traceID), self.getEarliest(traceID)],
                    [ax.get_ylim()[0],
                     ax.get_ylim()[1]],
                    'g:', label = 'earliest')
        if self.getLatest(traceID) is not None:
            ax.plot([self.getLatest(traceID), self.getLatest(traceID)],
                    [ax.get_ylim()[0],
                     ax.get_ylim()[1]],
                    'b:', label = 'latest')
        ax.plot([0, self.getPick(traceID)],
                 [folm * max(hoscf.getCF()), folm * max(hoscf.getCF())],
                 'm:', label = 'folm = %s' %folm)
        ax.set_xlabel('Time [s]')
        ax.legend()

    def plot3dttc(self, step = 0.5, contour = False, plotpicks = False, method = 'linear', ax = None):
        '''
        Plots a 3D 'traveltime cone' as surface plot by interpolating on a regular grid over the traveltimes, not yet regarding the vertical offset of the receivers.

        :param: step (optional), gives the stepsize for the interpolated grid. Default is 0.5
        :type: 'float'

        :param: contour (optional), plot contour plot instead of surface plot
        :type: 'logical'

        :param: plotpicks (optional), plot the data points onto the interpolated grid
        :type: 'logical'

        :param: method (optional), interpolation method; can be 'linear' (default) or 'cubic'
        :type: 'string'
        '''
        from scipy.interpolate import griddata
        from matplotlib import cm
        from mpl_toolkits.mplot3d import Axes3D

        x = []
        y = []
        z = []
        for traceID in self.picks.keys():
            if self.getPickFlag(traceID) != 0:
                x.append(self.getRecLoc(traceID)[0])
                y.append(self.getRecLoc(traceID)[1])
                z.append(self.getPick(traceID))

        xaxis = np.arange(min(x), max(x), step)
        yaxis = np.arange(min(y), max(y), step)
        xgrid, ygrid = np.meshgrid(xaxis, yaxis)
        zgrid = griddata((x, y), z, (xgrid, ygrid), method = method)

        if ax == None:
            fig = plt.figure()
            ax = plt.axes(projection = '3d')

        xsrc, ysrc, zsrc = self.getSrcLoc()

        if contour == True:
            ax.contour3D(xgrid,ygrid,zgrid,20)
        else:
            ax.plot_surface(xgrid, ygrid, zgrid, linewidth = 0, cmap = cm.jet, vmin = min(z), vmax = max(z))
        ax.plot([xsrc], [ysrc], [self.getPick(0)], 'k*', markersize = 20) # plot source location
        ax.plot([xsrc], [ysrc], [self.getPick(0)], 'r*', markersize = 15) # plot source location

        if plotpicks == True:
            ax.plot(x, y, z, 'k.')

    def plotttc(self, method, *args):
        plotmethod = {'2d': self.plot2dttc, '3d': self.plot3dttc}

        plotmethod[method](*args)

    def matshow(self, ax = None, step = 0.5, method = 'linear', plotRec = True, annotations = True, colorbar = True):
        '''
        Plots a 2D matrix of the interpolated traveltimes. This needs less performance than plot3dttc

        :param: step (optional), gives the stepsize for the interpolated grid. Default is 0.5
        :type: 'float'

        :param: method (optional), interpolation method; can be 'linear' (default) or 'cubic'
        :type: 'string'

        :param: plotRec (optional), plot the receiver positions (colored scatter plot, should not be
        deactivated because there might be receivers that are not inside the interpolated area)
        :type: 'logical'

        :param: annotations (optional), displays traceIDs as annotations
        :type: 'logical'
        '''
        from scipy.interpolate import griddata
        from matplotlib import cm
        cmap = cm.jet

        x = []; xcut = []
        y = []; ycut = []
        z = []; zcut = []

        for traceID in self.picks.keys():
            if self.getPickFlag(traceID) != 0:
                x.append(self.getRecLoc(traceID)[0])
                y.append(self.getRecLoc(traceID)[1])
                z.append(self.getPick(traceID))
            if self.getPickFlag(traceID) == 0 and self.getPickIncludeRemoved(traceID) is not None:
                xcut.append(self.getRecLoc(traceID)[0])
                ycut.append(self.getRecLoc(traceID)[1])
                zcut.append(self.getPickIncludeRemoved(traceID))

        tmin = 0.8 * min(z) # 20% cushion for colorbar
        tmax = 1.2 * max(z)

        xaxis = np.arange(min(x), max(x), step)
        yaxis = np.arange(min(y), max(y), step)
        xgrid, ygrid = np.meshgrid(xaxis, yaxis)
        zgrid = griddata((x, y), z, (xgrid, ygrid), method='linear')

        if ax == None:
            fig = plt.figure()
            ax = plt.axes()

        count = 0
        ax.imshow(zgrid, extent = [min(x), max(x), min(y), max(y)], vmin = tmin, vmax = tmax, cmap = cmap, origin = 'lower', alpha = 0.85)
        plt.text(0.45, 0.9, 'shot: %s' %self.getShotnumber(), transform = ax.transAxes)
        sc = ax.scatter(x, y, c = z, s = 30, label = 'picked shots', vmin = tmin, vmax = tmax, cmap = cmap, linewidths = 1.5)
        for xyz in zip(xcut, ycut, zcut):
            x, y, z = xyz
            label = None
            if z > tmax:
                count += 1
                z = 'w'
                if count == 1:
                    label = 'cut out shots'
            ax.scatter(x, y, c = z, s = 30, edgecolor = 'm', label = label, vmin = tmin, vmax = tmax, cmap = cmap, linewidths = 1.5)
        if colorbar == True:
            cbar = plt.colorbar(sc)
            cbar.set_label('Time [s]')

        ax.legend()
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.plot(self.getSrcLoc()[0], self.getSrcLoc()[1],'*k', markersize = 15) # plot source location

        if annotations == True:
            for traceID in self.getTraceIDlist():
                if self.getPickFlag(traceID) is not 0:
                    ax.annotate(' %s' %traceID , xy = (self.getRecLoc(traceID)[0], self.getRecLoc(traceID)[1]),
                                fontsize = 'x-small', color = 'k')
                else:
                    ax.annotate(' %s' %traceID , xy = (self.getRecLoc(traceID)[0], self.getRecLoc(traceID)[1]),
                                fontsize = 'x-small', color = 'r')

        plt.show()


