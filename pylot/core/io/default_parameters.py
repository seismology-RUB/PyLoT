#!/usr/bin/env python
# -*- coding: utf-8 -*-

defaults = {'rootpath': {'type': str,
                         'tooltip': 'project path',
                         'value': ''},
            
            'datapath': {'type': str,
                         'tooltip': 'data path',
                         'value': ''},
            
            'database': {'type': str,
                         'tooltip': 'name of data base',
                         'value': ''},
            
            'eventID': {'type': str,
                        'tooltip': 'event ID for single event processing (* for all events found in database)',
                        'value': ''},
            
            'extent': {'type': str,
                       'tooltip': 'extent of array ("local", "regional" or "global")',
                       'value': 'local'},
            
            'invdir': {'type': str,
                       'tooltip': 'full path to  inventory or dataless-seed file',
                       'value': ''},
            
            'datastructure': {'type': str,
                              'tooltip': 'choose data structure',
                              'value': 'PILOT'},
            
            'apverbose': {'type': bool,
                          'tooltip': "choose 'True' or 'False' for terminal output",
                          'value':  True},
            
            'nllocbin': {'type': str,
                         'tooltip': 'path to NLLoc executable',
                         'value': ''},
            
            'nllocroot': {'type': str,
                          'tooltip': 'root of NLLoc-processing directory',
                          'value': ''},
            
            'phasefile': {'type': str,
                          'tooltip': 'name of autoPyLoT-output phase file for NLLoc',
                          'value': 'AUTOPHASES.obs'},
            
            'ctrfile': {'type': str,
                        'tooltip': 'name of autoPyLoT-output control file for NLLoc',
                        'value': 'Insheim_min1d2015_auto.in'},
            
            'ttpatter': {'type': str,
                         'tooltip': 'pattern of NLLoc ttimes from grid',
                         'value': 'ttime'},
            
            'outpatter': {'type': str,
                          'tooltip': 'pattern of NLLoc-output file',
                          'value': 'AUTOLOC_nlloc'},
            
            'vp': {'type': float,
                   'tooltip': 'average P-wave velocity',
                   'value': 3530.},
            
            'rho': {'type': float,
                    'tooltip': 'average rock density [kg/m^3]',
                    'value': 2500.},
            
            'Qp': {'type': (float, float),
                   'tooltip': 'quality factor for P waves (Qp*f^a); list(Qp, a)',
                   'value': (300., 0.8)},
            
            'pstart': {'type': float,
                       'tooltip': 'start time [s] for calculating CF for P-picking',
                       'value': 15.0},
            
            'pstop': {'type': float,
                      'tooltip': 'end time [s] for calculating CF for P-picking',
                      'value': 60.0},
            
            'sstart': {'type': float,
                       'tooltip': 'start time [s] relative to P-onset for calculating CF for S-picking',
                       'value': -1.0},
            
            'sstop': {'type': float,
                      'tooltip': 'end time [s] after P-onset for calculating CF for S-picking',
                      'value': 10.0},
            
            'bpz1': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of first band pass filter Z-comp. [Hz]',
                     'value': (2, 20)},
            
            'bpz2': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of second band pass filter Z-comp. [Hz]',
                     'value': (2, 30)},
            
            'bph1': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of first band pass filter H-comp. [Hz]',
                     'value': (2, 15)},
            
            'bph2': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of second band pass filter z-comp. [Hz]',
                     'value': (2, 20)},
            
            'algoP': {'type': str,
                      'tooltip': 'choose algorithm for P-onset determination (HOS, ARZ, or AR3)',
                      'value': 'HOS'},
            
            'tlta': {'type': float,
                     'tooltip': 'for HOS-/AR-AIC-picker, length of LTA window [s]',
                     'value': 7.0},
            
            'hosorder': {'type': int,
                         'tooltip': 'for HOS-picker, order of Higher Order Statistics',
                         'value': 4},
            
            'Parorder': {'type': int,
                         'tooltip': 'for AR-picker, order of AR process of Z-component',
                         'value': 2},
            
            'tdet1z': {'type': float,
                       'tooltip': 'for AR-picker, length of AR determination window [s] for Z-component, 1st pick',
                       'value': 1.2},
            
            'tpred1z': {'type': float,
                        'tooltip': 'for AR-picker, length of AR prediction window [s] for Z-component, 1st pick',
                        'value': 0.4},
            
            'tdet2z': {'type': float,
                       'tooltip': 'for AR-picker, length of AR determination window [s] for Z-component, 2nd pick',
                       'value': 0.6},
            
            'tpred2z': {'type': float,
                        'tooltip': 'for AR-picker, length of AR prediction window [s] for Z-component, 2nd pick',
                        'value': 0.2},
            
            'addnoise': {'type': float,
                         'tooltip': 'add noise to seismogram for stable AR prediction',
                         'value': 0.001},
            
            'tsnrz': {'type': (float, float, float, float),
                      'tooltip': 'for HOS/AR, window lengths for SNR-and slope estimation [tnoise, tsafetey, tsignal, tslope] [s]',
                      'value': (3, 0.1, 0.5, 1.0)},
            
            'pickwinP': {'type': float,
                         'tooltip': 'for initial AIC pick, length of P-pick window [s]',
                         'value': 3.0},
            
            'Precalcwin': {'type': float,
                           'tooltip': 'for HOS/AR, window length [s] for recalculation of CF (relative to 1st pick)',
                           'value': 6.0},
            
            'aictsmooth': {'type': float,
                           'tooltip': 'for HOS/AR, take average of samples for smoothing of AIC-function [s]',
                           'value': 0.2},
            
            'tsmoothP': {'type': float,
                         'tooltip': 'for HOS/AR, take average of samples for smoothing CF [s]',
                         'value': 0.1},
            
            'ausP': {'type': float,
                     'tooltip': 'for HOS/AR, artificial uplift of samples (aus) of CF (P)',
                     'value': 0.001},
            
            'nfacP': {'type': float,
                      'tooltip': 'for HOS/AR, noise factor for noise level determination (P)',
                      'value': 1.3},
            
            'algoS': {'type': str,
                      'tooltip': 'choose algorithm for S-onset determination (ARH or AR3)',
                      'value': 'ARH'},
            
            'tdet1h': {'type': float,
                       'tooltip': 'for HOS/AR, length of AR-determination window [s], H-components, 1st pick',
                       'value': 0.8},
            
            'tpred1h': {'type': float,
                        'tooltip': 'for HOS/AR, length of AR-prediction window [s], H-components, 1st pick',
                        'value': 0.4},
            
            'tdet2h': {'type': float,
                       'tooltip': 'for HOS/AR, length of AR-determinaton window [s], H-components, 2nd pick',
                       'value': 0.6},
            
            'tpred2h': {'type': float,
                        'tooltip': 'for HOS/AR, length of AR-prediction window [s], H-components, 2nd pick',
                        'value': 0.3},
            
            'Sarorder': {'type': int,
                         'tooltip': 'for AR-picker, order of AR process of H-components',
                         'value': 4},
            
            'Srecalcwin': {'type': float,
                           'tooltip': 'for AR-picker, window length [s] for recalculation of CF (2nd pick) (H)',
                           'value': 5.0},
            
            'pickwinS': {'type': float,
                         'tooltip': 'for initial AIC pick, length of S-pick window [s]',
                         'value': 3.0},
            
            'tsnrh': {'type': (float, float, float, float),
                      'tooltip': 'for ARH/AR3, window lengths for SNR-and slope estimation [tnoise, tsafetey, tsignal, tslope] [s]',
                      'value': (2, 0.2, 1.5, 0.5)},
            
            'aictsmoothS': {'type': float,
                            'tooltip': 'for AIC-picker, take average of samples for smoothing of AIC-function [s]',
                            'value': 0.5},
            
            'tsmoothS': {'type': float,
                         'tooltip': 'for AR-picker, take average of samples for smoothing CF [s] (S)',
                         'value': 0.7},
            
            'ausS': {'type': float,
                     'tooltip': 'for HOS/AR, artificial uplift of samples (aus) of CF (S)',
                     'value': 0.9},
            
            'nfacS': {'type': float,
                      'tooltip': 'for AR-picker, noise factor for noise level determination (S)',
                      'value': 1.5},
            
            'minfmweight': {'type': int,
                            'tooltip': 'minimum required P weight for first-motion determination',
                            'value': 1},
            
            'minFMSNR': {'type': float,
                         'tooltip': 'miniumum required SNR for first-motion determination',
                         'value': 2.},
            
            'fmpickwin': {'type': float,
                          'tooltip': 'pick window around P onset for calculating zero crossings',
                          'value': 0.2},
            
            'timeerrorsP': {'type': (float, float, float, float),
                            'tooltip': 'discrete time errors [s] corresponding to picking weights [0 1 2 3] for P',
                            'value': (0.01, 0.02, 0.04, 0.08)},
            
            'timeerrorsS': {'type': (float, float, float, float),
                            'tooltip': 'discrete time errors [s] corresponding to picking weights [0 1 2 3] for S',
                            'value': (0.04, 0.08, 0.16, 0.32)},
            
            'minAICPslope': {'type': float,
                             'tooltip': 'below this slope [counts/s] the initial P pick is rejected',
                             'value': 0.8},
            
            'minAICPSNR': {'type': float,
                           'tooltip': 'below this SNR the initial P pick is rejected',
                           'value': 1.1},
            
            'minAICSslope': {'type': float,
                             'tooltip': 'below this slope [counts/s] the initial S pick is rejected',
                             'value': 1.},
            
            'minAICSSNR': {'type': float,
                           'tooltip': 'below this SNR the initial S pick is rejected',
                           'value': 1.5},
            
            'minsiglength': {'type': float,
                             'tooltip': 'length of signal part for which amplitudes must exceed noiselevel [s]',
                             'value': 1.},
            
            'noisefactor': {'type': float,
                            'tooltip': 'noiselevel*noisefactor=threshold',
                            'value': 1.0},
            
            'minpercent': {'type': float,
                           'tooltip': 'required percentage of amplitudes exceeding threshold',
                           'value': 10.},
            
            'zfac': {'type': float,
                     'tooltip': 'P-amplitude must exceed at least zfac times RMS-S amplitude',
                     'value': 1.5},
            
            'mdttolerance': {'type': float,
                             'tooltip': 'maximum allowed deviation of P picks from median [s]',
                             'value': 6.0},
            
            'wdttolerance': {'type': float,
                             'tooltip': 'maximum allowed deviation from Wadati-diagram',
                             'value': 1.0},
            
            'WAscaling': {'type': (float, float, float),
                         'tooltip': 'Scaling relation (log(Ao)+Alog(r)+Br+C) of Wood-Anderson amplitude Ao [nm] \
                                     If zeros are set, original Richter magnitude is calculated!',
                         'value': (0., 0., 0.)},

            'magscaling': {'type': (float, float),
                         'tooltip': 'Scaling relation for derived local magnitude [a*Ml+b]. \
                                     If zeros are set, no scaling of network magnitude is applied!',
                         'value': (0., 0.)}
}

settings_main={
    'dirs':[
        'rootpath',
        'datapath',
        'database',
        'eventID',
        'invdir',
        'datastructure',
        'apverbose'],
    'nlloc':[
        'nllocbin',
        'nllocroot',
        'phasefile',
        'ctrfile',
        'ttpatter',
        'outpatter'],
    'smoment':[
        'vp',
        'rho',
        'Qp'],
    'localmag':[
        'WAscaling',
        'magscaling'],
    'pick':[
        'extent',
        'pstart',
        'pstop',
        'sstart',
        'sstop',
        'bpz1',
        'bpz2',
        'bph1',
        'bph2']
}

settings_special_pick={
    'z':[
        'algoP',
        'tlta',
        'hosorder',
        'Parorder',
        'tdet1z',
        'tpred1z',
        'tdet2z',
        'tpred2z',
        'addnoise',
        'tsnrz',
        'pickwinP',
        'Precalcwin',
        'aictsmooth',
        'tsmoothP',
        'ausP',
        'nfacP'],
    'h':[
        'algoS',
        'tdet1h',
        'tpred1h',
        'tdet2h',
        'tpred2h',
        'Sarorder',
        'Srecalcwin',
        'pickwinS',
        'tsnrh',
        'aictsmoothS',
        'tsmoothS',
        'ausS',
        'nfacS'],
    'fm':[
        'minfmweight',
        'minFMSNR',
        'fmpickwin'],
    'quality':[
        'timeerrorsP',
        'timeerrorsS',
        'minAICPslope',
        'minAICPSNR',
        'minAICSslope',
        'minAICSSNR',
        'minsiglength',
        'noisefactor',
        'minpercent',
        'zfac',
        'mdttolerance',
        'wdttolerance']
}
