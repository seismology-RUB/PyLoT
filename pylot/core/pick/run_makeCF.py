#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
   Script to run autoPyLoT-script "makeCF.py".
   Only for test purposes!
"""

from obspy.core import read 
import matplotlib.pyplot as plt  
import numpy as np
from CharFuns import *
import glob
import argparse

def run_makeCF(project, database, event, iplot, station=None):
    #parameters for CF calculation
    t2 = 7              #length of moving window for HOS calculation [sec]
    p = 4               #order of statistics
    cuttimes = [10, 40] #start and end time vor CF calculation
    bpz = [2, 30]       #corner frequencies of bandpass filter, vertical component
    bph = [2, 15]       #corner frequencies of bandpass filter, horizontal components
    tdetz= 1.2          #length of AR-determination window [sec], vertical component
    tdeth= 0.8          #length of AR-determination window [sec], horizontal components
    tpredz = 0.4        #length of AR-prediction window [sec], vertical component
    tpredh = 0.4        #length of AR-prediction window [sec], horizontal components
    addnoise = 0.001    #add noise to seismogram for stable AR prediction
    arzorder = 2        #chosen order of AR process, vertical component
    arhorder = 4        #chosen order of AR process, horizontal components
    #get waveform data
    if station:
       dpz = '/DATA/%s/EVENT_DATA/LOCAL/%s/%s/%s*EHZ.msd' % (project, database, event, station)
       dpe = '/DATA/%s/EVENT_DATA/LOCAL/%s/%s/%s*EHE.msd' % (project, database, event, station)
       dpn = '/DATA/%s/EVENT_DATA/LOCAL/%s/%s/%s*EHN.msd' % (project, database, event, station)
    else:
       dpz = '/DATA/%s/EVENT_DATA/LOCAL/%s/%s/*EHZ.msd' % (project, database, event)
       dpe = '/DATA/%s/EVENT_DATA/LOCAL/%s/%s/*EHE.msd' % (project, database, event)
       dpn = '/DATA/%s/EVENT_DATA/LOCAL/%s/%s/*EHN.msd' % (project, database, event)
    wfzfiles = glob.glob(dpz)
    wfefiles = glob.glob(dpe)
    wfnfiles = glob.glob(dpn)
    if wfzfiles:
       for i in range(len(wfzfiles)):
          print 'Vertical component data found ...'
          print wfzfiles[i]
          st = read('%s' % wfzfiles[i])
          st_copy = st.copy()
          #filter and taper data
          tr_filt = st[0].copy()
          tr_filt.filter('bandpass', freqmin=bpz[0], freqmax=bpz[1], zerophase=False) 
          tr_filt.taper(max_percentage=0.05, type='hann')
          st_copy[0].data = tr_filt.data
          ##############################################################
          #calculate HOS-CF using subclass HOScf of class CharacteristicFunction
          hoscf = HOScf(st_copy, cuttimes, t2, p) #instance of HOScf
          ##############################################################
          #calculate AIC-HOS-CF using subclass AICcf of class CharacteristicFunction
          #class needs stream object => build it
          tr_aic = tr_filt.copy()
          tr_aic.data = hoscf.getCF()
          st_copy[0].data = tr_aic.data
          aiccf = AICcf(st_copy, cuttimes, t2, p) #instance of AICcf
          ##############################################################
          #calculate ARZ-CF using subclass ARZcf of class CharcteristicFunction
          #get stream object of filtered data
          st_copy[0].data = tr_filt.data
          arzcf = ARZcf(st_copy, cuttimes, tpredz, arzorder, tdetz, addnoise) #instance of ARZcf
          ##############################################################
          #calculate AIC-ARZ-CF using subclass AICcf of class CharacteristicFunction
          #class needs stream object => build it
          tr_arzaic = tr_filt.copy()
          tr_arzaic.data = arzcf.getCF()
          st_copy[0].data = tr_arzaic.data
          araiccf = AICcf(st_copy, cuttimes, t2, p) #instance of AICcf
    elif not wfzfiles:
       print 'No vertical component data found!'

    if wfefiles and wfnfiles:
       for i in range(len(wfefiles)):
          print 'Horizontal component data found ...'
          print wfefiles[i]
          print wfnfiles[i]
          #merge streams
          H = read('%s' % wfefiles[i])
          H += read('%s' % wfnfiles[i])
          H_copy = H.copy()
          #filter and taper data
          trH1_filt = H[0].copy()
          trH2_filt = H[1].copy()
          trH1_filt.filter('bandpass', freqmin=bph[0], freqmax=bph[1], zerophase=False) 
          trH2_filt.filter('bandpass', freqmin=bph[0], freqmax=bph[1], zerophase=False) 
          trH1_filt.taper(max_percentage=0.05, type='hann')
          trH2_filt.taper(max_percentage=0.05, type='hann')
          H_copy[0].data = trH1_filt.data
          H_copy[1].data = trH2_filt.data
          ##############################################################
          #calculate ARH-CF using subclass ARHcf of class CharcteristicFunction
          arhcf = ARHcf(H_copy, cuttimes, tpredh, arhorder, tdeth, addnoise) #instance of ARHcf
          ##############################################################
          if iplot:
             #plot vertical trace
             plt.figure()
             tr = st[0]
             #time vectors
             tdata = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
             tCF = np.arange(cuttimes[0], cuttimes[1], tr.stats.delta)
             tARZCF = np.arange(cuttimes[0] + tdetz + tpredz, cuttimes[1], tr.stats.delta)
             p1 = plt.plot(tdata, tr_filt.data/max(tr_filt.data), 'k')
             p2 = plt.plot(tCF, hoscf.getCF()/max(hoscf.getCF()), 'r')
             p3 = plt.plot(tCF, aiccf.getCF()/max(aiccf.getCF()), 'b')
             p4 = plt.plot(tARZCF, arzcf.getCF()/max(arzcf.getCF()), 'g')
             p5 = plt.plot(tARZCF, araiccf.getCF()/max(araiccf.getCF()), 'y')
             plt.yticks([])
             plt.xlabel('Time [s]')
             plt.ylabel('Normalized Counts')
             plt.title([tr.stats.station, tr.stats.channel])
             plt.suptitle(tr.stats.starttime)
             plt.legend([p1, p2, p3, p4, p5], ['Data', 'HOS-CF', 'HOSAIC-CF', 'ARZ-CF', 'ARZAIC-CF']) 
             #plot horizontal traces
             plt.figure(2)
             plt.subplot(211)
             th1data = np.arange(0, trH1_filt.stats.npts / trH1_filt.stats.sampling_rate, trH1_filt.stats.delta)
             th2data = np.arange(0, trH2_filt.stats.npts / trH2_filt.stats.sampling_rate, trH2_filt.stats.delta)
             tARHCF = np.arange(cuttimes[0] + tdeth + tpredh, cuttimes[1], trH1_filt.stats.delta)
             p21 = plt.plot(th1data, trH1_filt.data/max(trH1_filt.data), 'k')
             p22 = plt.plot(tARHCF, arhcf.getCF()/max(arhcf.getCF()), 'r') 
             plt.yticks([])
             plt.ylabel('Normalized Counts')
             plt.title([trH1_filt.stats.station, trH1_filt.stats.channel])
             plt.suptitle(trH1_filt.stats.starttime)
             plt.legend([p21, p22], ['Data', 'ARH-CF']) 
             plt.subplot(212)
             p23 = plt.plot(th2data, trH2_filt.data/max(trH2_filt.data), 'k')
             p24 = plt.plot(tARHCF, arhcf.getCF()/max(arhcf.getCF()), 'r') 
             plt.title([trH2_filt.stats.station, trH2_filt.stats.channel])
             plt.yticks([])
             plt.xlabel('Time [s]')
             plt.ylabel('Normalized Counts')
             plt.show()
             raw_input()
             plt.close()

parser = argparse.ArgumentParser()
parser.add_argument('--project', type=str, help='project name (e.g. Insheim)')
parser.add_argument('--database', type=str, help='event data base (e.g. 2014.09_Insheim)')
parser.add_argument('--event', type=str, help='event ID (e.g. e0010.015.14)')
parser.add_argument('--iplot', help='anything, if set, figure occurs')
parser.add_argument('--station', type=str, help='Station ID (e.g. INS3) (optional)')
args = parser.parse_args()

run_makeCF(args.project, args.database, args.event, args.iplot, args.station)
