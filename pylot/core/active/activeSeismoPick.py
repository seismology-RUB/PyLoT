import sys
import numpy as np
from pylot.core.active import seismicshot

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
        
    def _generateSurvey(self):
        from obspy.core import read

        shot_dict = {}
        shotlist = self.getShotlist()
        for shotnumber in shotlist:       # loop over data files
            # generate filenames and read manual picks to a list
            obsfile = self._obsdir + str(shotnumber) + '_pickle.dat'

            if not obsfile in shot_dict.keys():
                shot_dict[shotnumber] = []
            shot_dict[shotnumber] = seismicshot.SeismicShot(obsfile)
            shot_dict[shotnumber].setParameters('shotnumber', shotnumber)

        self.setArtificialPick(0, 0) # artificial pick at source origin

        self.data = shot_dict
        print ("Generated Survey object for %d shots" % len(shotlist))
        print ("Total number of traces: %d \n" %self.countAllTraces())

    def setParametersForShots(self, cutwindow = (0, 0.2), tmovwind = 0.3, tsignal = 0.03, tgap = 0.0007):
        if (cutwindow == (0, 0.2) and tmovwind == 0.3 and
            tsignal == 0.03 and tgap == 0.0007):
            print ("Warning: Standard values used for "
                   "setParamters. This may not be clever.")
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

                # ++ TEST: set and check SNR before adding to distance bin ############################
                shot.setSNR(traceID)    
                #if shot.getSNR(traceID)[0] < snrthreshold:
                if shot.getSNR(traceID)[0] < shot.getSNRthreshold(traceID):
                        shot.removePick(traceID)
                # -- TEST: set and check SNR before adding to distance bin ############################

                if shot.getPick(traceID) is not None:
                    shot.setEarllatepick(traceID)

            tpicksum += (datetime.now() - tstartpick); tpick = tpicksum/count
            tremain = (tpick * (len(self.getShotDict()) - count))
            tend = datetime.now() + tremain
            progress = float(count) / float(len(self.getShotDict())) * 100
            self._update_progress(shot.getShotname(), tend, progress)
        print('\npickAllShots: Finished\n')

    def setArtificialPick(self, traceID, pick):
        for shot in self.data.values():
            shot.setPick(traceID, pick)
            shot.setPickwindow(traceID, shot.getCut())

    def countAllTraces(self):
        numtraces = 0
        for line in self.getShotlist():
            for line in self.getReceiverlist():
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
                if shot.getPick(traceID) is not None:
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
                if shot.getPick(traceID) is not None:
                    pick = shot.getPick(traceID) * fmtomo_factor
                    delta = shot.getPickError(traceID) * fmtomo_factor
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
            if shot.getPick(traceID) is not None:
                count += 1
        return count

    def plotAllPicks(self, plotDeleted = False):
        '''
        Plots all picks over the distance between source and receiver. Returns (ax, region)
        '''
        import matplotlib.pyplot as plt
        import math
        plt.interactive(True)
        from pylot.core.active.surveyPlotTools import regions

        dist = []
        pick = []
        snrloglist = []
        for shot in self.data.values():
            for traceID in shot.getTraceIDlist():
                if plotDeleted == False:
                    if shot.getPick(traceID) is not None: 
                        dist.append(shot.getDistance(traceID))
                        pick.append(shot.getPick(traceID))
                        snrloglist.append(math.log10(shot.getSNR(traceID)[0]))
                elif plotDeleted == True:
                    dist.append(shot.getDistance(traceID))
                    pick.append(shot.getPick(traceID))
                    snrloglist.append(math.log10(shot.getSNR(traceID)[0]))

        ax = self.createPlot(dist, pick, snrloglist, label = 'log10(SNR)')
        region = regions(ax, self.data)
        ax.legend()

        return ax, region

    def createPlot(self, dist, pick, inkByVal, label):
        import matplotlib.pyplot as plt
        plt.interactive(True)
        cm = plt.cm.jet

        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig = ax.scatter(dist, pick, cmap = cm, c = inkByVal, s = 5, edgecolors = 'none', label = label)
        cbar = plt.colorbar(fig, fraction = 0.05)
        cbar.set_label(label)
        plt.title('Plot of all Picks')
        plt.xlabel('Distance [m]')
        plt.ylabel('Time [s]')

        return ax

    def _update_progress(self, shotname, tend, progress):
        sys.stdout.write("Working on shot %s. ETC is %02d:%02d:%02d [%2.2f %%]\r" 
                         %(shotname, tend.hour, tend.minute, tend.second, progress))
        sys.stdout.flush()

    def saveSurvey(self, filename = 'survey.pickle'):
        import cPickle
        outfile = open(filename, 'wb')

        cPickle.dump(self, outfile, -1)
        print('saved Survey to file %s'%(filename))
        
    @staticmethod    
    def from_pickle(filename):
        import cPickle
        infile = open(filename, 'rb')
        return cPickle.load(infile)
