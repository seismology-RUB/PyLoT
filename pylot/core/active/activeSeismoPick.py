from pylot.core.active import seismicshot
import numpy as np

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

    def pickAllShots(self, HosAic = 'hos', vmin = 333, vmax = 5500):
        '''
        Automatically pick all traces of all shots of the survey.
        '''
        from datetime import datetime
        count = 0

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
                shot.pickTraces(traceID, pickmethod, windowsize, folm, HosAic) # picker

                # ++ TEST: set and check SNR before adding to distance bin ############################
                shot.setSNR(traceID)    
                #if shot.getSNR(traceID)[0] < snrthreshold:
                if shot.getSNR(traceID)[0] < shot.getSNRthreshold(traceID):
                        shot.removePick(traceID)
                # -- TEST: set and check SNR before adding to distance bin ############################

                if shot.getPick(traceID) is not None:
                    shot.setEarllatepick(traceID)

            tpicksum += (datetime.now() - tstartpick); tpick = tpicksum/count
            tremain = (tpick * (len(survey.getShotDict()) - count))
            tend = datetime.now() + tremain
            print 'shot: %s, est. time to be finished is %s:%s:%s' % (shot.getShotname(), tend.hour, tend.minute, tend.second)




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
