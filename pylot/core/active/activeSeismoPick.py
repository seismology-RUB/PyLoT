# -*- coding: utf-8 -*-
import sys
import numpy as np
from pylot.core.active import seismicshot
from pylot.core.active.surveyUtils import cleanUp

class Survey(object):
    def __init__(self, path, sourcefile, receiverfile, useDefaultParas = False):
        '''
        The Survey Class contains all shots [type: seismicshot] of a survey
        as well as the aquisition geometry and the topography.
        '''
        self.data = {}
        self._topography = None
        self._recfile = receiverfile
        self._sourcefile = sourcefile
        self._obsdir = path
        self._generateSurvey()
        if useDefaultParas == True:
            self.setParametersForShots()
        self._removeAllEmptyTraces()
        self._updateShots()
        self.setArtificialPick(0, 0)

    def _generateSurvey(self):
        from obspy.core import read

        shot_dict = {}
        shotlist = self.getShotlist()
        for shotnumber in shotlist:       # loop over data files
            # generate filenames and read manual picks to a list
            obsfile = self._obsdir + str(shotnumber) + '_pickle.dat'
            if obsfile not in shot_dict.keys():
                shot_dict[shotnumber] = []
            shot_dict[shotnumber] = seismicshot.SeismicShot(obsfile)
            shot_dict[shotnumber].setParameters('shotnumber', shotnumber)

        self.data = shot_dict
        print ("Generated Survey object for %d shots" % len(shotlist))
        print ("Total number of traces: %d \n" %self.countAllTraces())

    def setArtificialPick(self, traceID, pick):
        '''
        Sets an artificial pick for a traceID of all shots in the survey object.
        (This can be used to create a pick with t = 0 at the source origin)
        '''
        for shot in self.data.values():
            shot.setPick(traceID, pick)

    def setParametersForShots(self, cutwindow = (0, 0.2), tmovwind = 0.3, tsignal = 0.03, tgap = 0.0007):
        if (cutwindow == (0, 0.2) and tmovwind == 0.3 and
            tsignal == 0.03 and tgap == 0.0007):
            print ("Warning: Standard values used for "
                   "setParamters. This might not be clever.")
        # CHANGE this later. Parameters only needed for survey, not for each shot.
        for shot in self.data.values():
            shot.setCut(cutwindow)
            shot.setTmovwind(tmovwind)
            shot.setTsignal(tsignal)
            shot.setTgap(tgap)
            shot.setRecfile(self.getPath() + self.getReceiverfile())
            shot.setSourcefile(self.getPath() + self.getSourcefile())
            shot.setOrder(order = 4)
        print ("setParametersForShots: Parameters set to:\n"
               "cutwindow = %s, tMovingWindow = %f, tsignal = %f, tgap = %f"
               %(cutwindow, tmovwind, tsignal, tgap))

    def _removeAllEmptyTraces(self):
        filename = 'removeEmptyTraces.out'
        count = 0
        for shot in self.data.values():
            removed = shot.removeEmptyTraces()
            if removed is not None:
                if count == 0: outfile = open(filename, 'w')
                count += 1
                outfile.writelines('shot: %s, removed empty traces: %s\n'
                                   %(shot.getShotnumber(), removed))
        print ("\nremoveEmptyTraces: Finished! Removed %d traces" %count)
        if count > 0:
            print ("See %s for more information "
                   "on removed traces."%(filename))
            outfile.close()

    def _updateShots(self):
        filename = 'updateShots.out'
        count = 0; countTraces = 0
        for shot in self.data.values():
            del_traceIDs = shot.updateTraceList()
            if len(del_traceIDs) > 0:
                if count == 0: outfile = open(filename, 'w')
                count += 1
                countTraces += len(del_traceIDs)
                outfile.writelines("shot: %s, removed traceID(s) %s because "
                                   "they were not found in the corresponding stream\n"
                                   %(shot.getShotnumber(), del_traceIDs))

        print ("\nupdateShots: Finished! Updated %d shots and removed "
               "%d traces" %(count, countTraces))
        if count > 0:
            print ("See %s for more information "
                   "on removed traces."%(filename))
            outfile.close()

    def pickAllShots(self, windowsize, HosAic = 'hos', vmin = 333, vmax = 5500, folm = 0.6):
        '''
        Automatically pick all traces of all shots of the survey.
        '''
        from datetime import datetime
        starttime = datetime.now()
        count = 0; tpicksum = starttime - starttime

        for shot in self.data.values():
            tstartpick = datetime.now(); count += 1
            for traceID in shot.getTraceIDlist():
                distance = shot.getDistance(traceID) # receive distance

                pickwin_used = shot.getCut()
                cutwindow = shot.getCut()

                # for higher distances use a linear vmin/vmax to cut out late/early regions with high noise
                if distance > 5.:
                    pwleft = distance/vmax ################## TEST
                    pwright = distance/vmin
                    if pwright > cutwindow[1]:
                        pwright = cutwindow[1]
                    pickwin_used = (pwleft, pwright)

                shot.setPickwindow(traceID, pickwin_used)
                shot.pickTraces(traceID, windowsize, folm, HosAic) # picker

                shot.setSNR(traceID)
                #if shot.getSNR(traceID)[0] < snrthreshold:
                if shot.getSNR(traceID)[0] < shot.getSNRthreshold(traceID):
                        shot.removePick(traceID)

                # set epp and lpp if SNR > 1 (else earllatepicker cant set values)
                if shot.getSNR(traceID)[0] > 1:
                    shot.setEarllatepick(traceID)

            tpicksum += (datetime.now() - tstartpick); tpick = tpicksum/count
            tremain = (tpick * (len(self.getShotDict()) - count))
            tend = datetime.now() + tremain
            progress = float(count) / float(len(self.getShotDict())) * 100
            self._update_progress(shot.getShotname(), tend, progress)
        print('\npickAllShots: Finished\n')
        ntraces = self.countAllTraces()
        pickedtraces = self.countAllPickedTraces()
        print('Picked %s / %s traces (%d %%)\n'
              %(pickedtraces, ntraces, float(pickedtraces)/float(ntraces)*100.))


    def recover(self):
        '''
        Recovers all (accidently) removed picks. Still regards SNR threshold.
        '''
        print('Recovering survey...')
        numpicks = 0
        for shot in self.data.values():
            for traceID in shot.getTraceIDlist():
                if shot.getFlag(traceID) == 0:
                    shot.setFlag(traceID, 1)
                    if shot.getSNR(traceID)[0] < shot.getSNRthreshold(traceID):
                        shot.removePick(traceID)
                    else:
                        numpicks += 1
        print('Recovered %d picks'%numpicks)

    def setArtificialPick(self, traceID, pick):
        for shot in self.data.values():
            shot.setPick(traceID, pick)
            shot.setPickwindow(traceID, shot.getCut())

    def countAllTraces(self):
        numtraces = 0
        for shot in self.getShotlist():
            for rec in self.getReceiverlist(): ### shot.getReceiverlist etc.
                numtraces += 1
        return numtraces

    def getShotlist(self):
        filename = self.getPath() + self.getSourcefile()
        srcfile = open(filename,'r')
        shotlist = []
        for line in srcfile.readlines():
            line = line.split()
            shotlist.append(int(line[0]))

        return shotlist

    def getReceiverlist(self):
        filename = self.getPath() + self.getReceiverfile()
        recfile = open(filename,'r')
        reclist = []
        for line in recfile.readlines():
            line = line.split()
            reclist.append(int(line[0]))

        return reclist

    def getShotDict(self):
        return self.data

    def getShot(self, shotnumber):
        return self.data[shotnumber]

    def getSourcefile(self):
        return self._sourcefile

    def getReceiverfile(self):
        return self._recfile

    def getPath(self):
        return self._obsdir

    def getStats(self):
        info_dict = {}
        for shot in self.data.values():
            pickedTraces = 0
            snrlist = []
            dist = []
            numtraces = len(shot.getTraceIDlist())
            for traceID in shot.getTraceIDlist():
                snrlist.append(shot.getSNR(traceID)[0])
                dist.append(shot.getDistance(traceID))
                if shot.getFlag(traceID) is not 0:
                    pickedTraces += 1
            info_dict[shot.getShotnumber()] = {'numtraces': numtraces,
                                               'picked traces': [pickedTraces,
                                                                 '%2.2f %%'%(float(pickedTraces) /
                                                                             float(numtraces) * 100)],
                                               'mean SNR': np.mean(snrlist),
                                               'mean distance': np.mean(dist)}

        return info_dict

    def getShotForShotnumber(self, shotnumber):
        for shot in self.data.values():
            if shot.getShotnumber() == shotnumber:
                return shot

    def exportFMTOMO(self, directory = 'FMTOMO_export', sourcefile = 'input_sf.in', ttFileExtension = '.tt'):
        def getAngle(distance):
            PI = np.pi
            R = 6371.
            angle = distance * 180 / (PI * R)
            return angle

        count = 0
        fmtomo_factor = 1000 # transforming [m/s] -> [km/s]
        LatAll = []; LonAll = []; DepthAll = []
        srcfile = open(directory + '/' + sourcefile, 'w')
        srcfile.writelines('%10s\n' %len(self.data)) # number of sources
        for shotnumber in self.getShotlist():
            shot = self.getShotForShotnumber(shotnumber)
            ttfilename = str(shotnumber) + ttFileExtension
            (x, y, z) = shot.getSrcLoc() # getSrcLoc returns (x, y, z)
            srcfile.writelines('%10s %10s %10s\n' %(getAngle(y), getAngle(x), (-1)*z)) # lat, lon, depth
            LatAll.append(getAngle(y)); LonAll.append(getAngle(x)); DepthAll.append((-1)*z)
            srcfile.writelines('%10s\n' %1) #
            srcfile.writelines('%10s %10s %10s\n' %(1, 1, ttfilename))
            ttfile = open(directory + '/' + ttfilename, 'w')
            traceIDlist = shot.getTraceIDlist()
            traceIDlist.sort()
            ttfile.writelines(str(self.countPickedTraces(shot)) + '\n')
            for traceID in traceIDlist:
                if shot.getFlag(traceID) is not 0:
                    pick = shot.getPick(traceID) * fmtomo_factor
                    delta = shot.getSymmetricPickError(traceID) * fmtomo_factor
                    (x, y, z) = shot.getRecLoc(traceID)
                    ttfile.writelines('%20s %20s %20s %10s %10s\n' %(getAngle(y), getAngle(x), (-1)*z, pick, delta))
                    LatAll.append(getAngle(y)); LonAll.append(getAngle(x)); DepthAll.append((-1)*z)
                    count += 1
            ttfile.close()
        srcfile.close()
        print 'Wrote output for %s traces' %count
        print 'WARNING: output generated for FMTOMO-obsdata. Obsdata seems to take Lat, Lon, Depth and creates output for FMTOMO as Depth, Lat, Lon'
        print 'Dimensions of the seismic Array, transformed for FMTOMO, are Depth(%s, %s), Lat(%s, %s), Lon(%s, %s)'%(
            min(DepthAll), max(DepthAll), min(LatAll), max(LatAll), min(LonAll), max(LonAll))

    def countPickedTraces(self, shot):
        count = 0
        for traceID in shot.getTraceIDlist():
            if shot.getFlag(traceID) is not 0:
                count += 1
        return count

    def countAllPickedTraces(self):
        count = 0
        for shot in self.data.values():
            for traceID in shot.getTraceIDlist():
                if shot.getFlag(traceID) is not 0:
                    count += 1
        return count

    def plotAllShots(self, rows = 3, columns = 4):
        '''
        Plots all shots as Matrices with the color corresponding to the traveltime for each receiver.
        IMPORTANT NOTE: Topography (z - coordinate) is not considered in the diagrams!
        '''
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        plt.interactive(True)

        fig = plt.figure()
        ax = fig.add_subplot(111)

        figPerSubplot = columns * rows

        index = 1
        #shotnames = []
        #shotnumbers = []

        # for shot in self.data.values():
        #     shotnames.append(shot.getShotname())
        #     shotnumbers.append(shot.getShotnumber())

        # shotnumbers = [shotnumbers for (shotnumbers, shotnames) in sorted(zip(shotnumbers, shotnames))]

        for shotnumber in self.getShotlist():
            if index <= figPerSubplot:
                #ax = fig.add_subplot(3,3,i, projection = '3d', title = 'shot:'
                #+str(shot_dict[shotnumber].getShotnumber()), xlabel = 'X', ylabel = 'Y', zlabel = 'traveltime')
                #shot_dict[shotnumber].plot3dttc(ax = ax, plotpicks = True)
                ax = fig.add_subplot(3, 4, index)
                self.getShot(shotnumber).matshow(ax = ax, colorbar = False, annotations = True)
                index += 1
            if index > figPerSubplot:
                fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, wspace = 0, hspace = 0)
                fig = plt.figure()
                index = 1

        fig.subplots_adjust(left = 0, bottom = 0, right = 1, top = 1, wspace = 0, hspace = 0)

    def plotAllPicks(self, plotRemoved = False, colorByVal = 'log10SNR', ax = None, cbar = None, refreshPlot = False):
        '''
        Plots all picks over the distance between source and receiver. Returns (ax, region).
        Picks can be checked and removed by using region class (pylot.core.active.surveyPlotTools.regions)

        :param: plotRemoved, if True plots traces that were picked but removed from the survey (flag = 0)
        :type: logical

        :param: colorByVal, can be "log10SNR", "pickerror", or "spe"
        :type: str

        Examples:

        regions.chooseRectangles():
         - lets the user choose several rectangular regions in the plot

        regions.plotTracesInRegions():
         - creates plots (shot.plot_traces) for all traces in the active regions (i.e. chosen by e.g. chooseRectangles)

        regions.setActiveRegionsForDeletion():
         - highlights all shots in a the active regions for deletion

        regions.deleteMarkedPicks():
         - deletes the picks (pick flag set to 0) for all shots set for deletion

        regions.deselectSelection(number):
         - deselects the region of number = number

        '''

        import matplotlib.pyplot as plt
        import math
        plt.interactive(True)
        from pylot.core.active.surveyPlotTools import regions

        dist = []
        pick = []
        snrlog = []
        pickerror = []
        spe = []

        for shot in self.data.values():
            for traceID in shot.getTraceIDlist():
                if plotRemoved == False:
                    if shot.getFlag(traceID) is not 0 or plotRemoved == True:
                        dist.append(shot.getDistance(traceID))
                        pick.append(shot.getPick(traceID))
                        snrlog.append(math.log10(shot.getSNR(traceID)[0]))
                        pickerror.append(shot.getPickError(traceID))
                        spe.append(shot.getSymmetricPickError(traceID))

        color = {'log10SNR': snrlog,
         'pickerror': pickerror,
         'spe': spe}
        self.color = color
        if refreshPlot is False:
            ax, cbar = self.createPlot(dist, pick, color[colorByVal], label='%s' % colorByVal)
            region = regions(ax, cbar, self)
            ax.legend()
            return (ax, region)
        if refreshPlot is True:
            ax, cbar = self.createPlot(dist, pick, color[colorByVal], label='%s' % colorByVal, ax=ax, cbar=cbar)
            ax.legend()
            return ax

    def createPlot(self, dist, pick, inkByVal, label, ax = None, cbar = None):
        import matplotlib.pyplot as plt
        plt.interactive(True)
        cm = plt.cm.jet
        if ax is None:
            print('Generating new plot...')
            fig = plt.figure()
            ax = fig.add_subplot(111)
            sc = ax.scatter(dist, pick, cmap=cm, c=inkByVal, s=5, edgecolors='none', label=label)
            cbar = plt.colorbar(sc, fraction=0.05)
            cbar.set_label(label)
            ax.set_xlabel('Distance [m]')
            ax.set_ylabel('Time [s]')
            ax.text(0.5, 0.95, 'Plot of all picks', transform=ax.transAxes, horizontalalignment='center')
        else:
            sc = ax.scatter(dist, pick, cmap=cm, c=inkByVal, s=5, edgecolors='none', label=label)
            cbar = plt.colorbar(sc, cax=cbar.ax)
            cbar.set_label(label)
            ax.set_xlabel('Distance [m]')
            ax.set_ylabel('Time [s]')
            ax.text(0.5, 0.95, 'Plot of all picks', transform=ax.transAxes, horizontalalignment='center')
        return (ax, cbar)

    def _update_progress(self, shotname, tend, progress):
        sys.stdout.write('Working on shot %s. ETC is %02d:%02d:%02d [%2.2f %%]\r' % (shotname,
         tend.hour,
         tend.minute,
         tend.second,
         progress))
        sys.stdout.flush()

    def saveSurvey(self, filename = 'survey.pickle'):
        import cPickle
        cleanUp(self)
        outfile = open(filename, 'wb')

        cPickle.dump(self, outfile, -1)
        print('saved Survey to file %s'%(filename))

    @staticmethod
    def from_pickle(filename):
        import cPickle
        infile = open(filename, 'rb')
        return cPickle.load(infile)
