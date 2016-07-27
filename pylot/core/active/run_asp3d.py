#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
from obspy import read
from obspy import Stream
from obspy import Trace
from datetime import datetime
import numpy as np

from pylot.core.active import surveyUtils
from pylot.core.active import seismicshot
from pylot.core.active import activeSeismoPick
from pylot.core.active import fmtomoUtils
from pylot.core.active import seismicArrayPreparation
# reload(seismicshot)
# reload(surveyUtils)
# reload(activeSeismoPick)

#####################################################################################
# parameter definitions:#############################################################
cutwindow = (0, 0.15)              # cut out a part of the trace [seconds]
tmovwind = 0.1                     # size of the moving window
windowsize = (30, 0)               # windowsize for AIC picker (indices around HOS picks [time/sampling rate] !!!)
folm = 0.6                         # fraction of local maximum for threshold picker
tsignal = 0.03
tgap = 0.0007

nproc = 1

vmin = 333
vmax = 5500

HosAic = 'aic' # can be 'hos' or 'aic'

rockeskyll = False
GZB = True
bb1 = False

# Simulation parameters #############################################################
simulation = True

niter = 4

bottomBoundary = -50.0
topBoundary = 5.0
nPointsPropgrid = (100, 100, 100)
nPointsVgrid = (30, 30, 30)
cushionfactor = 0.1
interpolationMethod = 'linear'
mygrid = '../mygrid.in'
                   
cwd = os.getcwd()
simuldir = 'fmtomo_simulation'
pickdir = 'picks'
fmtomopath = '/rscratch/minos22/marcel/flachseismik/fmtomo/GZB_clean/'
######################################################################################
######################################################################################
if rockeskyll == True:
    receiverfile = "Geophone_interpoliert_rockes"
    sourcefile = "Schusspunkte_rockes"
    obsdir = "/rscratch/minos22/marcel/flachseismik/rockeskyll_200615_270615/"
    filename = 'survey_rockes.pickle'
elif GZB == True:
    receiverfile = "Geophone_interpoliert_GZB"
    sourcefile = "Schusspunkte_GZB"
    #sourcefile = "Schusspunkte_GZB_1shot"
    obsdir = "/rscratch/minos22/marcel/flachseismik/GZB_26_06_15_01/"
    filename = 'survey_GZB.pickle'
elif bb1 == True:
    receiverfile = "Geophone_Marcel" 
    sourcefile = "Schusspunkte_Marcel"
    obsdir = "/rscratch/minos22/marcel/flachseismik/bachelor_bausenberg/"
    filename = 'survey_bb1.pickle'
######################################################################################

starttime = datetime.now()

print('\n--------------------Starting Script at %s -------------------\n' %starttime.time())
print('directory: %s\nsourcefile: %s\nreceiverfile: %s\nsurvey output filename: %s\n' %(obsdir, sourcefile, receiverfile, filename))
if HosAic == 'aic': print('picking with AIC\n')
if HosAic == 'hos': print('picking with HOS\n')

survey = activeSeismoPick.Survey(obsdir, os.path.join(obsdir, sourcefile),
                                 os.path.join(obsdir, receiverfile), useDefaultParas = False)
survey.setParametersForAllShots(cutwindow, tmovwind, tsignal, tgap)
surveyUtils.setDynamicFittedSNR(survey.getShotDict())
#surveyUtils.setConstantSNR(survey.getShotDict(), 0)
survey.setArtificialPick(0, 0) # artificial pick at source origin
print('\nDone after %s seconds!\n------------------------------------------------------------------------------\n' % (datetime.now() - starttime).seconds)
survey.pickAllShots(vmin, vmax, folm, HosAic, windowsize, cores = nproc)
survey.cleanBySPE(maxSPE = 0.0075)
survey.saveSurvey(filename)
print('\n--- Finished picking ---')
print'Elapsed time:', datetime.now()-starttime


######################################################################################
if simulation == False:
    sys.exit()

survey = activeSeismoPick.Survey.from_pickle(filename)

if not os.path.isdir(os.path.join(cwd, simuldir)):
    err = os.mkdir(os.path.join(cwd, simuldir))

picks = os.path.join(simuldir, pickdir)

if not os.path.isdir(os.path.join(cwd, picks)):
    err2 = os.mkdir(os.path.join(cwd, picks))

survey.exportFMTOMO(picks)

####### hard coded
os.chdir(simuldir)
survey.loadArray(obsdir, 'Geophone_eingemessen_GZB', 'Schusspunkte_GZB')
survey.seisArray.generateFMTOMOinputFromArray(nPointsPropgrid, nPointsVgrid,
                                              (bottomBoundary, topBoundary), cushionfactor,
                                              interpolationMethod, customgrid = mygrid,
                                              writeVTK = False)

os.chdir(cwd)
####### test

tomo = fmtomoUtils.Tomo3d(fmtomopath, simuldir)
tomo.runTOMO3D(nproc, niter)

print('\n--- Finished script ---')
print'Elapsed time:', datetime.now()-starttime
