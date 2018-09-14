# AutoPyLoT Tuning

A description of the parameters used for determing automatic picks.

## Filter parameters and cut times

Parameters applied to the traces before picking algorithm starts.

Name | Description
--- | ---
*P Start*, *P Stop* | Define time interval relative to trace start time for CF calculation on vertical trace. Value is relative to theoretical onset time if 'Use TayPy' option is enabled in main settings of 'Tune Autopicker' dialog.
*S Start*, *S Stop* | Define time interval relative to trace start time for CF calculation on horizontal traces. Value is relative to theoretical onset time if 'Use TayPy' option is enabled in main settings of 'Tune Autopicker' dialog.
*Bandpass Z1* | Filter settings for Butterworth bandpass applied to vertical trace for calculation of initial P pick.
*Bandpass Z2* | Filter settings for Butterworth bandpass applied to vertical trace for calculation of precise P pick.
*Bandpass H1* | Filter settings for Butterworth bandpass applied to horizontal traces for calculation of initial S pick.
*Bandpass H2* | Filter settings for Butterworth bandpass applied to horizontal traces for calculation of precise S pick.

## Inital P pick

Parameters used for determination of initial P pick.

Name | Description
--- | ---
*tLTA* | Size of gliding LTA window in seconds used for calculation of HOS-CF.
*pickwin P* | Size of time window in seconds in which the minimum of the AIC-CF in front of the maximum of the HOS-CF is determined.
*AICtsmooth* | Average of samples in this time window will be used for smoothing of the AIC-CF.
*checkwinP* | Time in front of the global maximum of the HOS-CF in which to search for a second local extrema.
*minfactorP* | Used with *checkwinP*. If a second local maximum is found, it has to be at least as big as the first maximum * *minfactorP*.
*tsignal* | Time window in seconds after the initial P pick used for determining signal amplitude.
*tnoise* | Time window in seconds in front of initial P pick used for determining noise amplitude.
*tsafetey* | Time in seconds between *tsignal* and *tnoise*.
*tslope* | Time window in seconds after initial P pick in which the slope of the onset is calculated.

## Inital S pick

Parameters used for determination of initial S pick

Name | Description
--- | ---
*tdet1h* | Length of time window in seconds in which AR params of the waveform are determined.
*tpred1h* | Length of time window in seconds in which the waveform is predicted using the AR modell.
*AICtsmoothS* | Average of samples in this time window is used for smoothing the AIC-CF.
*pickwinS* | Time window in which the minimum in the AIC-CF in front of the maximum in the ARH cf is determined.
*checkwinS* | Time in front of the global maximum of the ARH-CF in which to search for a second local extrema.
*minfactorP* | Used with *checkwinS*. If a second local maximum is found, it has to be at least as big as the first maximum * *minfactorS*.
*tsignal* | Time window in seconds after the initial P pick used for determining signal amplitude.
*tnoise* | Time window in seconds in front of initial P pick used for determining noise amplitude.
*tsafetey* | Time in seconds between *tsignal* and *tnoise*.
*tslope* | Time window in seconds after initial P pick in which the slope of the onset is calculated.

## Precise P pick

Parameters used for determination of precise P pick.

Name | Description
--- | ---
*Precalcwin* | Time window in seconds for recalculation of the HOS-CF. The new CF will be two times the size of *Precalcwin*, since it will be calculated from the intialpick to +/- *Precalcwin*.
*tsmoothP* | Average of samples in this time window will be used for smoothing the secod HOS-CF.
*ausP* | Controls artificial uplift of samples during precise picking. A common local minimum of the smoothed and unsmoothed HOS-CF is found when the previous sample is larger or equal to the current sample times (1+*ausP*).

## Precise S pick

Parameters used for determination of precise S pick.

Name | Description
--- | ---
*tdet2h* | Time window for determination of AR coefficients.
*tpred2h* | Time window in which the waveform is predicted using the determined AR parameters.
*Srecalcwin* | Time window for recalculation of ARH-CF. New CF will be calculated from initial pick +/- *Srecalcwin*.
*tsmoothS* | Average of samples in this time window will be used for smoothing the second ARH-CF.
*ausS* | Controls artificial uplift of samples during precise picking. A common local minimum of the smoothed and unsmoothed ARH-CF is found when the previous sample is larger or equal to the current sample times (1+*ausS*). 
*pickwinS* | Time window around intial pick in which to look for a precise pick.


## Pick quality control

Parameters used for checking quality and integrity of automatic picks.

Name | Description
--- | ---
*minAICPslope* | Initial P picks with a slope lower than this value will be discared.
*minAICPSNR* | Initial P picks with a SNR below this value will be discarded.
*minAICSslope* | Initial S picks with a slope lower than this value will be discarded.
*minAICSSNR* | Initial S picks with a SNR below this value will be discarded.
*minsiglength*, *noisefacor*. *minpercent* | Parameters for checking signal length. In the time window of size *minsiglength* after the initial P pick *minpercent* of samples have to be larger than the RMS value.
*zfac* | To recognize misattributed S picks, the RMS amplitude of vertical and horizontal traces are compared. The RMS amplitude of the vertical traces has to be at least *zfac* higher than the RMS amplitude on the horizontal traces for the pick to be accepted as a valid P pick.
*jackfactor* | A P pick is removed if the jackknife pseudo value of the variance of his subgroup ist larger than the variance of all picks multiplied with the *jackfactor*.
*mdttolerance* | Maximum allowed deviation of P onset times from the median .
*wdttolerance* | Maximum allowed deviation of S onset times from the line during the Wadati test.

## Pick quality determination

Parameters for discrete quality classes.

Name | Description
--- | ---
*timeerrorsP* | Width of the time windows in seconds between earliest and latest possible pick which represent the quality classes 0, 1, 2, 3 for P onsets.
*timeerrorsS* |  Width of the time windows in seconds between earliest and latest possible pick which represent the quality classes 0, 1, 2, 3 for S onsets.
*nfacP*, *nfacS* | For determination of latest possible onset time. The time when the signal reaches an amplitude of *nfac* * mean value of the RMS amplitude in the time window *tnoise* corresponds to the latest possible onset time.

