import numpy as np

def readParameters(parfile, parameter):
    from ConfigParser import ConfigParser
    parameterConfig = ConfigParser()
    parameterConfig.read('parfile')

    value = parameterConfig.get('vars', parameter).split('\t')[0]

    return value

def setArtificialPick(shot_dict, traceID, pick):
    for shot in shot_dict.values():
        shot.setPick(traceID, pick)
        shot.setPickwindow(traceID, shot.getCut())

def fitSNR4dist(shot_dict, shiftdist = 30, shiftSNR = 100):
    import numpy as np
    import matplotlib.pyplot as plt
    dists = []
    picks = []
    snrs = []
    snr_sqrt_inv = []
    snrthresholds = []
    snrBestFit = []
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist():
            if shot.getSNR(traceID)[0] >= 1:
                dists.append(shot.getDistance(traceID))
                picks.append(shot.getPickIncludeRemoved(traceID))
                snrs.append(shot.getSNR(traceID)[0])
                snr_sqrt_inv.append(1/np.sqrt(shot.getSNR(traceID)[0]))
    fit = np.polyfit(dists, snr_sqrt_inv, 1)
    fit_fn = np.poly1d(fit)
    for dist in dists:
        snrBestFit.append((1/(fit_fn(dist)**2)))
        dist += shiftdist
        snrthresholds.append((1/(fit_fn(dist)**2)) - shiftSNR * np.exp(-0.05 * dist))
    plotFittedSNR(dists, snrthresholds, snrs, snrBestFit)
    return fit_fn #### ZU VERBESSERN, sollte fertige funktion wiedergeben


def plotFittedSNR(dists, snrthresholds, snrs, snrBestFit):
    import matplotlib.pyplot as plt
    plt.interactive(True)
    fig = plt.figure()
    plt.plot(dists, snrs, 'b.', markersize = 2.0, label = 'SNR values')
    dists.sort()
    snrthresholds.sort(reverse = True)
    snrBestFit.sort(reverse = True)
    plt.plot(dists, snrthresholds, 'r', markersize = 1, label = 'Fitted threshold')
    plt.plot(dists, snrBestFit, 'k', markersize = 1, label = 'Best fitted curve')
    plt.xlabel('Distance[m]')
    plt.ylabel('SNR')
    plt.legend()

def setDynamicFittedSNR(shot_dict, shiftdist = 30, shiftSNR = 100, p1 = 0.004, p2 = -0.0007):
    import numpy as np
    minSNR = 2.5
    #fit_fn = fitSNR4dist(shot_dict)
    fit_fn = np.poly1d([p1, p2])
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist(): ### IMPROVE
            dist = shot.getDistance(traceID) + shiftdist
            snrthreshold = (1/(fit_fn(dist)**2)) - shiftSNR * np.exp(-0.05 * dist)
            if snrthreshold < minSNR:
                print('WARNING: SNR threshold %s lower %s. Set SNR threshold to %s.'
                      %(snrthreshold, minSNR, minSNR))
                shot.setSNRthreshold(traceID, minSNR)
            else:
                shot.setSNRthreshold(traceID, snrthreshold)
    print "setDynamicFittedSNR: Finished setting of fitted SNR-threshold"

def setConstantSNR(shot_dict, snrthreshold = 2.5):
    import numpy as np
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist():
            shot.setSNRthreshold(traceID, snrthreshold)
    print "setConstantSNR: Finished setting of SNR threshold to a constant value of %s"%snrthreshold

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

def cleanUp(survey):

    for shot in survey.data.values():
        shot.traces4plot = {}
