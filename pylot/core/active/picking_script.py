import sys
from obspy import read
from obspy import Stream
from obspy import Trace
from datetime import datetime
import numpy as np

from pylot.core.active import surveyUtils
from pylot.core.active import seismicshot
import activeSeismoPick
reload(seismicshot)
reload(surveyUtils)

#####################################################################################
# parameter definitions:#############################################################
#traceslist = list(np.arange(1, 49, 1))      # traces (1-48)
#shotlist = list(np.arange(302, 352, 1))   # shot-files(.dat) (Paffrath: 302-352) (Hauburg: 353-401) (arange+1!)
cutwindow = (0, 0.2)              # cut out a part of the trace [seconds]
tmovwind = 0.3                     # size of the moving window
windowsize = (5, 0)               # windowsize for AIC picker (indices around HOS picks!!)
pickwindow = cutwindow              # for local max and threshold picker: fraction of the seismogram used (0...1) TO BE DONE: depends on cutwindow!!!!
folm = 0.6

rockeskyll = False
if rockeskyll == True:
    receiverfile = "Geophone_interpoliert_rockes"
    sourcefile = "Schusspunkte_rockes"
    obsdir = "../rockeskyll_200615_270615/"
else:
    receiverfile = "Geophone_interpoliert_GZB" 
    sourcefile = "Schusspunkte_GZB"
    obsdir = "../GZB_26_06_15_01/"

# SNR
tsignal = 0.03
tgap = 0.0007
snrthreshold = 1
######################################################################################

vmin = 333
vmax = 5500
distBinSize = 2

###########################################
############# Settings: ###################
thresholdpick=True
localmaxpick=False

if thresholdpick == True: pickmethod = "threshold"
if localmaxpick == True: pickmethod = "localmax"

HosAic = 'hos' # can be 'hos' or 'aic'
###########################################

starttime = datetime.now()

print '\n--------------------Starting Script at %s -------------------\n' %starttime.time()
if thresholdpick == True: print 'Using treshold pickmethod!\n'
elif localmaxpick == True: print 'Using local maximum pick method!\n'
if HosAic == 'aic': print 'picking with AIC\n'
if HosAic == 'hos': print 'picking with HOS\n'

survey = activeSeismoPick.Survey(obsdir, sourcefile, receiverfile, True)
surveyUtils.setFittedSNR(survey.getShotDict())
print '\nDone after %s seconds!\n------------------------------------------------------------------------------\n' % (datetime.now() - starttime).seconds

count = 0; tpicksum = starttime - starttime

for shot in survey.data.values():
    tstartpick = datetime.now(); count += 1
    for traceID in shot.getTraceIDlist():
        distance = shot.getDistance(traceID) # receive distance

        pickwin_used = pickwindow # use pickwindow set in the parameter section
        # for higher distances use a linear vmin/vmax to cut out late/early regions with high noise
        if distance > 5.:
            pwleft = distance/vmax ################## TEST
            pwright = distance/vmin
            if pwright > cutwindow[1]:
                pwright = cutwindow[1]
        pickwin_used = (pwleft, pwright)

        shot.setPickwindow(traceID, pickwin_used)
        shot.pickTraces(traceID, windowsize, folm, HosAic) # picker
        #shot.setManualPicks(traceID, picklist) # set manual picks if given (yet used on 2D only)
        
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

print '\n--- Finished script ---'
print 'Elapsed time:', datetime.now()-starttime
