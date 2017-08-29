#!/usr/bin/env python
# -*- coding: utf-8 -*-

defaults = {'rootpath': {'type': str,
                         'tooltip': 'project path',
                         'value': '',
                         'namestring': 'Root path'},

            'datapath': {'type': str,
                         'tooltip': 'data path',
                         'value': '',
                         'namestring': 'Data path'},

            'database': {'type': str,
                         'tooltip': 'name of data base',
                         'value': '',
                         'namestring': 'Database path'},

            'eventID': {'type': str,
                        'tooltip': 'event ID for single event processing (* for all events found in database)',
                        'value': '',
                        'namestring': 'Event ID'},

            'extent': {'type': str,
                       'tooltip': 'extent of array ("local", "regional" or "global")',
                       'value': 'local',
                       'namestring': 'Array extent'},

            'invdir': {'type': str,
                       'tooltip': 'full path to  inventory or dataless-seed file',
                       'value': '',
                       'namestring': 'Inversion dir'},

            'datastructure': {'type': str,
                              'tooltip': 'choose data structure',
                              'value': 'PILOT',
                              'namestring': 'Datastructure'},

            'apverbose': {'type': bool,
                          'tooltip': "choose 'True' or 'False' for terminal output",
                          'value': True,
                          'namestring': 'App. verbosity'},

            'nllocbin': {'type': str,
                         'tooltip': 'path to NLLoc executable',
                         'value': '',
                         'namestring': 'NLLoc bin path'},

            'nllocroot': {'type': str,
                          'tooltip': 'root of NLLoc-processing directory',
                          'value': '',
                          'namestring': 'NLLoc root path'},

            'phasefile': {'type': str,
                          'tooltip': 'name of autoPyLoT-output phase file for NLLoc',
                          'value': 'AUTOPHASES.obs',
                          'namestring': 'Phase filename'},

            'ctrfile': {'type': str,
                        'tooltip': 'name of autoPyLoT-output control file for NLLoc',
                        'value': 'Insheim_min1d2015_auto.in',
                        'namestring': 'Control filename'},

            'ttpatter': {'type': str,
                         'tooltip': 'pattern of NLLoc ttimes from grid',
                         'value': 'ttime',
                         'namestring': 'Traveltime pattern'},

            'outpatter': {'type': str,
                          'tooltip': 'pattern of NLLoc-output file',
                          'value': 'AUTOLOC_nlloc',
                          'namestring': 'NLLoc output pattern'},

            'vp': {'type': float,
                   'tooltip': 'average P-wave velocity',
                   'value': 3530.,
                   'namestring': 'P-velocity'},

            'rho': {'type': float,
                    'tooltip': 'average rock density [kg/m^3]',
                    'value': 2500.,
                    'namestring': 'Density'},

            'Qp': {'type': (float, float),
                   'tooltip': 'quality factor for P waves (Qp*f^a); list(Qp, a)',
                   'value': (300., 0.8),
                   'namestring': ('Quality factor', 'Qp1', 'Qp2')},

            'pstart': {'type': float,
                       'tooltip': 'start time [s] for calculating CF for P-picking (if TauPy:'
                                  ' seconds relative to estimated onset)',
                       'value': 15.0,
                       'namestring': 'P start'},

            'pstop': {'type': float,
                      'tooltip': 'end time [s] for calculating CF for P-picking (if TauPy:'
                                 ' seconds relative to estimated onset)',
                      'value': 60.0,
                      'namestring': 'P stop'},

            'sstart': {'type': float,
                       'tooltip': 'start time [s] relative to P-onset for calculating CF for S-picking',
                       'value': -1.0,
                       'namestring': 'S start'},

            'sstop': {'type': float,
                      'tooltip': 'end time [s] after P-onset for calculating CF for S-picking',
                      'value': 10.0,
                      'namestring': 'S stop'},

            'bpz1': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of first band pass filter Z-comp. [Hz]',
                     'value': (2, 20),
                     'namestring': ('Z-bandpass 1', 'Lower', 'Upper')},

            'bpz2': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of second band pass filter Z-comp. [Hz]',
                     'value': (2, 30),
                     'namestring': ('Z-bandpass 2', 'Lower', 'Upper')},

            'bph1': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of first band pass filter H-comp. [Hz]',
                     'value': (2, 15),
                     'namestring': ('H-bandpass 1', 'Lower', 'Upper')},

            'bph2': {'type': (float, float),
                     'tooltip': 'lower/upper corner freq. of second band pass filter z-comp. [Hz]',
                     'value': (2, 20),
                     'namestring': ('H-bandpass 2', 'Lower', 'Upper')},

            'algoP': {'type': str,
                      'tooltip': 'choose algorithm for P-onset determination (HOS, ARZ, or AR3)',
                      'value': 'HOS',
                      'namestring': 'P algorithm'},

            'tlta': {'type': float,
                     'tooltip': 'for HOS-/AR-AIC-picker, length of LTA window [s]',
                     'value': 7.0,
                     'namestring': 'LTA window'},

            'hosorder': {'type': int,
                         'tooltip': 'for HOS-picker, order of Higher Order Statistics',
                         'value': 4,
                         'namestring': 'HOS order'},

            'Parorder': {'type': int,
                         'tooltip': 'for AR-picker, order of AR process of Z-component',
                         'value': 2,
                         'namestring': 'AR order P'},

            'tdet1z': {'type': float,
                       'tooltip': 'for AR-picker, length of AR determination window [s] for Z-component, 1st pick',
                       'value': 1.2,
                       'namestring': 'AR det. window Z 1'},

            'tpred1z': {'type': float,
                        'tooltip': 'for AR-picker, length of AR prediction window [s] for Z-component, 1st pick',
                        'value': 0.4,
                        'namestring': 'AR pred. window Z 1'},

            'tdet2z': {'type': float,
                       'tooltip': 'for AR-picker, length of AR determination window [s] for Z-component, 2nd pick',
                       'value': 0.6,
                       'namestring': 'AR det. window Z 2'},

            'tpred2z': {'type': float,
                        'tooltip': 'for AR-picker, length of AR prediction window [s] for Z-component, 2nd pick',
                        'value': 0.2,
                        'namestring': 'AR pred. window Z 2'},

            'addnoise': {'type': float,
                         'tooltip': 'add noise to seismogram for stable AR prediction',
                         'value': 0.001,
                         'namestring': 'Add noise'},

            'tsnrz': {'type': (float, float, float, float),
                      'tooltip': 'for HOS/AR, window lengths for SNR-and slope estimation [tnoise, tsafetey, tsignal, tslope] [s]',
                      'value': (3, 0.1, 0.5, 1.0),
                      'namestring': ('SNR windows P', 'Noise', 'Safety', 'Signal', 'Slope')},

            'pickwinP': {'type': float,
                         'tooltip': 'for initial AIC pick, length of P-pick window [s]',
                         'value': 3.0,
                         'namestring': 'AIC window P'},

            'Precalcwin': {'type': float,
                           'tooltip': 'for HOS/AR, window length [s] for recalculation of CF (relative to 1st pick)',
                           'value': 6.0,
                           'namestring': 'Recal. window P'},

            'aictsmooth': {'type': float,
                           'tooltip': 'for HOS/AR, take average of samples for smoothing of AIC-function [s]',
                           'value': 0.2,
                           'namestring': 'AIC smooth P'},

            'tsmoothP': {'type': float,
                         'tooltip': 'for HOS/AR, take average of samples for smoothing CF [s]',
                         'value': 0.1,
                         'namestring': 'CF smooth P'},

            'ausP': {'type': float,
                     'tooltip': 'for HOS/AR, artificial uplift of samples (aus) of CF (P)',
                     'value': 0.001,
                     'namestring': 'Artificial uplift P'},

            'nfacP': {'type': float,
                      'tooltip': 'for HOS/AR, noise factor for noise level determination (P)',
                      'value': 1.3,
                      'namestring': 'Noise factor P'},

            'algoS': {'type': str,
                      'tooltip': 'choose algorithm for S-onset determination (ARH or AR3)',
                      'value': 'ARH',
                      'namestring': 'S algorithm'},

            'tdet1h': {'type': float,
                       'tooltip': 'for HOS/AR, length of AR-determination window [s], H-components, 1st pick',
                       'value': 0.8,
                       'namestring': 'AR det. window H 1'},

            'tpred1h': {'type': float,
                        'tooltip': 'for HOS/AR, length of AR-prediction window [s], H-components, 1st pick',
                        'value': 0.4,
                        'namestring': 'AR pred. window H 1'},

            'tdet2h': {'type': float,
                       'tooltip': 'for HOS/AR, length of AR-determinaton window [s], H-components, 2nd pick',
                       'value': 0.6,
                       'namestring': 'AR det. window H 2'},

            'tpred2h': {'type': float,
                        'tooltip': 'for HOS/AR, length of AR-prediction window [s], H-components, 2nd pick',
                        'value': 0.3,
                        'namestring': 'AR pred. window H 2'},

            'Sarorder': {'type': int,
                         'tooltip': 'for AR-picker, order of AR process of H-components',
                         'value': 4,
                         'namestring': 'AR order S'},

            'Srecalcwin': {'type': float,
                           'tooltip': 'for AR-picker, window length [s] for recalculation of CF (2nd pick) (H)',
                           'value': 5.0,
                           'namestring': 'Recal. window S'},

            'pickwinS': {'type': float,
                         'tooltip': 'for initial AIC pick, length of S-pick window [s]',
                         'value': 3.0,
                         'namestring': 'AIC window S'},

            'tsnrh': {'type': (float, float, float, float),
                      'tooltip': 'for ARH/AR3, window lengths for SNR-and slope estimation [tnoise, tsafetey, tsignal, tslope] [s]',
                      'value': (2, 0.2, 1.5, 0.5),
                      'namestring': ('SNR windows S', 'Noise', 'Safety', 'Signal', 'Slope')},

            'aictsmoothS': {'type': float,
                            'tooltip': 'for AIC-picker, take average of samples for smoothing of AIC-function [s]',
                            'value': 0.5,
                            'namestring': 'AIC smooth S'},

            'tsmoothS': {'type': float,
                         'tooltip': 'for AR-picker, take average of samples for smoothing CF [s] (S)',
                         'value': 0.7,
                         'namestring': 'CF smooth S'},

            'ausS': {'type': float,
                     'tooltip': 'for HOS/AR, artificial uplift of samples (aus) of CF (S)',
                     'value': 0.9,
                     'namestring': 'Artificial uplift S'},

            'nfacS': {'type': float,
                      'tooltip': 'for AR-picker, noise factor for noise level determination (S)',
                      'value': 1.5,
                      'namestring': 'Noise factor S'},

            'minfmweight': {'type': int,
                            'tooltip': 'minimum required P weight for first-motion determination',
                            'value': 1,
                            'namestring': 'Min. P weight'},

            'minFMSNR': {'type': float,
                         'tooltip': 'miniumum required SNR for first-motion determination',
                         'value': 2.,
                         'namestring': 'Min SNR'},

            'fmpickwin': {'type': float,
                          'tooltip': 'pick window around P onset for calculating zero crossings',
                          'value': 0.2,
                          'namestring': 'Zero crossings window'},

            'timeerrorsP': {'type': (float, float, float, float),
                            'tooltip': 'discrete time errors [s] corresponding to picking weights [0 1 2 3] for P',
                            'value': (0.01, 0.02, 0.04, 0.08),
                            'namestring': ('Time errors P', '0', '1', '2', '3')},

            'timeerrorsS': {'type': (float, float, float, float),
                            'tooltip': 'discrete time errors [s] corresponding to picking weights [0 1 2 3] for S',
                            'value': (0.04, 0.08, 0.16, 0.32),
                            'namestring': ('Time errors S', '0', '1', '2', '3')},

            'minAICPslope': {'type': float,
                             'tooltip': 'below this slope [counts/s] the initial P pick is rejected',
                             'value': 0.8,
                             'namestring': 'Min. slope P'},

            'minAICPSNR': {'type': float,
                           'tooltip': 'below this SNR the initial P pick is rejected',
                           'value': 1.1,
                           'namestring': 'Min. SNR P'},

            'minAICSslope': {'type': float,
                             'tooltip': 'below this slope [counts/s] the initial S pick is rejected',
                             'value': 1.,
                             'namestring': 'Min. slope S'},

            'minAICSSNR': {'type': float,
                           'tooltip': 'below this SNR the initial S pick is rejected',
                           'value': 1.5,
                           'namestring': 'Min. SNR S'},

            'minsiglength': {'type': float,
                             'tooltip': 'length of signal part for which amplitudes must exceed noiselevel [s]',
                             'value': 1.,
                             'namestring': 'Min. signal length'},

            'noisefactor': {'type': float,
                            'tooltip': 'noiselevel*noisefactor=threshold',
                            'value': 1.0,
                            'namestring': 'Noise factor'},

            'minpercent': {'type': float,
                           'tooltip': 'required percentage of amplitudes exceeding threshold',
                           'value': 10.,
                           'namestring': 'Min amplitude [%]'},

            'zfac': {'type': float,
                     'tooltip': 'P-amplitude must exceed at least zfac times RMS-S amplitude',
                     'value': 1.5,
                     'namestring': 'Z factor'},

            'mdttolerance': {'type': float,
                             'tooltip': 'maximum allowed deviation of P picks from median [s]',
                             'value': 6.0,
                             'namestring': 'Median tolerance'},

            'wdttolerance': {'type': float,
                             'tooltip': 'maximum allowed deviation from Wadati-diagram',
                             'value': 1.0,
                             'namestring': 'Wadati tolerance'},

            'jackfactor': {'type': float,
                             'tooltip': 'pick is removed if the variance of the subgroup with the pick removed is larger than the mean variance of all subgroups times safety factor',
                             'value': 5.0,
                             'namestring': 'Jackknife safety factor'},

            'WAscaling': {'type': (float, float, float),
                          'tooltip': 'Scaling relation (log(Ao)+Alog(r)+Br+C) of Wood-Anderson amplitude Ao [nm] \
                          If zeros are set, original Richter magnitude is calculated!',
                          'value': (0., 0., 0.),
                          'namestring': ('Wood-Anderson scaling', '', '', '')},

            'magscaling': {'type': (float, float),
                           'tooltip': 'Scaling relation for derived local magnitude [a*Ml+b]. \
                           If zeros are set, no scaling of network magnitude is applied!',
                           'value': (0., 0.),
                           'namestring': ('Local mag. scaling', '', '')},

            'minfreq': {'type': (float, float),
                        'tooltip': 'Lower filter frequency [P, S]',
                        'value': (1.0, 1.0),
                        'namestring': ('Lower freq.', 'P', 'S')},

            'maxfreq': {'type': (float, float),
                        'tooltip': 'Upper filter frequency [P, S]',
                        'value': (10.0, 10.0),
                        'namestring': ('Upper freq.', 'P', 'S')},

            'filter_order': {'type': (int, int),
                             'tooltip': 'filter order [P, S]',
                             'value': (2, 2),
                             'namestring': ('Order', 'P', 'S')},

            'filter_type': {'type': (str, str),
                            'tooltip': 'filter type (bandpass, bandstop, lowpass, highpass) [P, S]',
                            'value': ('bandpass', 'bandpass'),
                            'namestring': ('Type', 'P', 'S')},

            'use_taup': {'type': bool,
                         'tooltip': 'use estimated traveltimes from TauPy for calculating windows for CF',
                         'value': True,
                         'namestring': 'Use TauPy'},

            'taup_model': {'type': str,
                           'tooltip': 'define TauPy model for traveltime estimation. Possible values: 1066a, 1066b, ak135, ak135f, herrin, iasp91, jb, prem, pwdk, sp6',
                           'value': 'iasp91',
                           'namestring': 'TauPy model'}
            }

settings_main = {
    'dirs': [
        'rootpath',
        'datapath',
        'database',
        'eventID',
        'invdir',
        'datastructure',
        'apverbose'],
    'nlloc': [
        'nllocbin',
        'nllocroot',
        'phasefile',
        'ctrfile',
        'ttpatter',
        'outpatter'],
    'smoment': [
        'vp',
        'rho',
        'Qp'],
    'localmag': [
        'WAscaling',
        'magscaling'],
    'filter': [
        'minfreq',
        'maxfreq',
        'filter_order',
        'filter_type'],
    'pick': [
        'extent',
        'pstart',
        'pstop',
        'sstart',
        'sstop',
        'use_taup',
        'taup_model',
        'bpz1',
        'bpz2',
        'bph1',
        'bph2']
}

settings_special_pick = {
    'z': [
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
    'h': [
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
    'fm': [
        'minfmweight',
        'minFMSNR',
        'fmpickwin'],
    'quality': [
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
        'wdttolerance',
        'jackfactor'],
}
