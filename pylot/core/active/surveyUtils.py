import numpy as np

def readParameters(parfile, parameter):
    from ConfigParser import ConfigParser
    parameterConfig = ConfigParser()
    parameterConfig.read('parfile')
    
    value = parameterConfig.get('vars', parameter).split('#')[0]
    value = value.replace(" ", "")

    return value

def fitSNR4dist(shot_dict, shiftdist = 5):
    dists = []
    picks = []
    snrs = []
    snr_sqrt_inv = []
    snrthresholds = []
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist():
            if shot.getSNR(traceID)[0] >= 1:
                dists.append(shot.getDistance(traceID))
                picks.append(shot.getPick_backup(traceID))
                snrs.append(shot.getSNR(traceID)[0])
                snr_sqrt_inv.append(1/np.sqrt(shot.getSNR(traceID)[0]))
    fit = np.polyfit(dists, snr_sqrt_inv, 1)
    fit_fn = np.poly1d(fit)
    for dist in dists:
        dist += shiftdist
        snrthresholds.append(1/(fit_fn(dist)**2))
    plotFittedSNR(dists, snrthresholds, snrs)
    return fit_fn #### ZU VERBESSERN, sollte fertige funktion wiedergeben

def plotFittedSNR(dists, snrthresholds, snrs):
    import matplotlib.pyplot as plt
    plt.interactive(True)
    fig = plt.figure()
    plt.plot(dists, snrs, '.', markersize = 0.5, label = 'SNR values')
    plt.plot(dists, snrthresholds, 'r.', markersize = 1, label = 'Fitted threshold')
    plt.xlabel('Distance[m]')
    plt.ylabel('SNR')
    plt.legend()

def setFittedSNR(shot_dict, shiftdist = 5, p1 = 0.004, p2 = -0.004):
    #fit_fn = fitSNR4dist(shot_dict)
    fit_fn = np.poly1d([p1, p2])
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist(): ### IMPROVE
            shot.setSNRthreshold(traceID, 1/(fit_fn(shot.getDistance(traceID) + shiftdist)**2)) ### s.o.
    print "\nsetFittedSNR: Finished setting of fitted SNR-threshold"

def findTracesInRanges(shot_dict, distancebin, pickbin):
    '''
    Returns traces corresponding to a certain area in a plot with all picks over the distances.
    
    :param: shot_dict, dictionary containing all shots that are used
    :type: dictionary
                                                                                     
    :param: distancebin
    :type: tuple, (dist1[m], dist2[m])
                                                                                     
    :param: pickbin
    :type: tuple, (t1[s], t2[s])
    '''
    shots_found = {}
    for shot in shot_dict.values():
        if shot.getTraceIDs4Dist(distancebin = distancebin) is not None:
            for traceID in shot.getTraceIDs4Dist(distancebin = distancebin):
                if pickbin[0] < shot.getPick(traceID) < pickbin[1]:
                    if shot.getShotnumber() not in shots_found.keys():
                        shots_found[shot.getShotnumber()] = []
                    shots_found[shot.getShotnumber()].append(traceID)

    return shots_found
