from __future__ import print_function


def readParameters(parfile, parameter):
    """

    :param parfile:
    :param parameter:
    :return:
    """
    from ConfigParser import ConfigParser
    parameterConfig = ConfigParser()
    parameterConfig.read('parfile')

    value = parameterConfig.get('vars', parameter).split('\t')[0]

    return value


def fitSNR4dist(shot_dict, shiftdist=30, shiftSNR=100):
    """
    Approach to fit the decreasing SNR with wave travel distance.

    :param shot_dict: dictionary containing Seismicshot objects (e.g. survey.getShotDict())
    :param shiftdist: shift compensating curve by a certain distance to the left
    :param shiftSNR: shift compensating curve by a certain SNR value to the bottom
    :return:
    """
    import numpy as np
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
                snr_sqrt_inv.append(1 / np.sqrt(shot.getSNR(traceID)[0]))
    fit = np.polyfit(dists, snr_sqrt_inv, 1)
    fit_fn = np.poly1d(fit)
    for dist in dists:
        snrBestFit.append((1 / (fit_fn(dist) ** 2)))
        dist += shiftdist
        snrthresholds.append((1 / (fit_fn(dist) ** 2)) - shiftSNR * np.exp(-0.05 * dist))
    plotFittedSNR(dists, snrthresholds, snrs, snrBestFit)
    return fit_fn  #### ZU VERBESSERN, sollte fertige funktion wiedergeben


def plotFittedSNR(dists, snrthresholds, snrs, snrBestFit):
    """

    :param dists:
    :param snrthresholds:
    :param snrs:
    :param snrBestFit:
    :return:
    """
    import matplotlib.pyplot as plt
    plt.interactive(True)
    fig = plt.figure()
    plt.plot(dists, snrs, 'b.', markersize=2.0, label='SNR values')
    dists.sort()
    snrthresholds.sort(reverse=True)
    snrBestFit.sort(reverse=True)
    plt.plot(dists, snrthresholds, 'r', markersize=1, label='Fitted threshold')
    plt.plot(dists, snrBestFit, 'k', markersize=1, label='Best fitted curve')
    plt.xlabel('Distance[m]')
    plt.ylabel('SNR')
    plt.legend()


def setDynamicFittedSNR(shot_dict, shiftdist=30, shiftSNR=100, p1=0.004, p2=-0.0007):
    """
    Set SNR values for a dictionary containing Seismicshots (e.g. survey.getShotDict())
    by parameters calulated from fitSNR4dist.

    :param shot_dict:
    :type shot_dict: dict
    :param shiftdist:
    :type shiftdist: int
    :param shiftSNR:
    :type shiftSNR: int
    :param p1:
    :type p1: float
    :param p2:
    :type p2: float
    :return:
    """
    import numpy as np
    minSNR = 2.5
    # fit_fn = fitSNR4dist(shot_dict)
    fit_fn = np.poly1d([p1, p2])
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist():  ### IMPROVE
            dist = shot.getDistance(traceID) + shiftdist
            snrthreshold = (1 / (fit_fn(dist) ** 2)) - shiftSNR * np.exp(-0.05 * dist)
            if snrthreshold < minSNR:
                print('WARNING: SNR threshold %s lower %s. Set SNR threshold to %s.'
                      % (snrthreshold, minSNR, minSNR))
                shot.setSNRthreshold(traceID, minSNR)
            else:
                shot.setSNRthreshold(traceID, snrthreshold)
    print("setDynamicFittedSNR: Finished setting of fitted SNR-threshold")


def setConstantSNR(shot_dict, snrthreshold=2.5):
    """
    Set a constant SNR value to all Seismicshots in a dictionary (e.g. survey.getShotDict()).

    :param shot_dict:
    :param snrthreshold:
    :return:
    """
    for shot in shot_dict.values():
        for traceID in shot.getTraceIDlist():
            shot.setSNRthreshold(traceID, snrthreshold)
    print("setConstantSNR: Finished setting of SNR threshold to a constant value of %s" % snrthreshold)


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
        if shot.getTraceIDs4Dist(distancebin=distancebin) is not None:
            for traceID in shot.getTraceIDs4Dist(distancebin=distancebin):
                if pickbin[0] < shot.getPick(traceID) < pickbin[1]:
                    if shot.getShotnumber() not in shots_found.keys():
                        shots_found[shot.getShotnumber()] = []
                    shots_found[shot.getShotnumber()].append(traceID)

    return shots_found


def cleanUp(survey):
    """
    Cleans up a Survey object by removing frontend information
    which can not be saved in the pickle format.
    """
    for shot in survey.data.values():
        shot.traces4plot = {}

def plotScatterStats4Shots(survey, variable, ax = None, twoDim = False):
    """
    Statistics, scatter plot.

    :param survey:
    :param variable: can be 'mean SNR', 'median SNR', 'mean SPE', 'median SPE', or 'picked traces'
    :return:
    """
    import matplotlib.pyplot as plt
    import numpy as np
    statsShot = {}
    x = []
    y = []
    value = []
    for shot in survey.data.values():
        for traceID in shot.getTraceIDlist():
            if not shot in statsShot.keys():
                statsShot[shot] = {'x': shot.getSrcLoc()[0],
                                   'y': shot.getSrcLoc()[1],
                                   'SNR': [],
                                   'SPE': [],
                                   'picked traces': 0}

            SNR = shot.getSNR(traceID)[0]
            if not SNR == np.inf:
                statsShot[shot]['SNR'].append(SNR)
            if shot.getPickFlag(traceID) == 1:
                statsShot[shot]['picked traces'] += 1
                statsShot[shot]['SPE'].append(shot.getSymmetricPickError(traceID))

    for shot in statsShot.keys():
        statsShot[shot]['mean SNR'] = np.mean(statsShot[shot]['SNR'])
        statsShot[shot]['median SNR'] = np.median(statsShot[shot]['SNR'])
        statsShot[shot]['mean SPE'] = np.mean(statsShot[shot]['SPE'])
        statsShot[shot]['median SPE'] = np.median(statsShot[shot]['SPE'])

    for shot in statsShot.keys():
        x.append(statsShot[shot]['x'])
        y.append(statsShot[shot]['y'])
        value.append(statsShot[shot][variable])

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    size = []
    for val in value:
        size.append(100 * val / max(value))

    sc = ax.scatter(x, y, s=size, c=value)
    ax.text(0.5, 1.05,'Plot of all shots',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    if not twoDim:
        ax.set_aspect('equal')
    cbar = ax.figure.colorbar(sc)
    cbar.set_label(variable)

    for shot in statsShot.keys():
        ax.annotate(' %s' % shot.getShotnumber(), xy=(shot.getSrcLoc()[0], shot.getSrcLoc()[1]),
                    fontsize='x-small', color='k')


def plotScatterStats4Receivers(survey, variable, ax = None, twoDim = False):
    """
    Statistics, scatter plot.

    :param survey:
    :param variable: can be 'mean SNR', 'median SNR', 'mean SPE', 'median SPE', or 'picked traces'
    :return:
    """
    import matplotlib.pyplot as plt
    import numpy as np
    statsRec = {}
    x = []
    y = []
    value = []
    for shot in survey.data.values():
        for traceID in shot.getTraceIDlist():
            if not traceID in statsRec.keys():
                statsRec[traceID] = {'x': shot.getRecLoc(traceID)[0],
                                     'y': shot.getRecLoc(traceID)[1],
                                     'SNR': [],
                                     'SPE': [],
                                     'picked traces': 0}

            SNR = shot.getSNR(traceID)[0]
            if not SNR == np.inf:
                statsRec[traceID]['SNR'].append(SNR)
            if shot.getPickFlag(traceID) == 1:
                statsRec[traceID]['picked traces'] += 1
                statsRec[traceID]['SPE'].append(shot.getSymmetricPickError(traceID))

    for traceID in statsRec.keys():
        statsRec[traceID]['mean SNR'] = np.mean(statsRec[traceID]['SNR'])
        statsRec[traceID]['median SNR'] = np.median(statsRec[traceID]['SNR'])
        statsRec[traceID]['mean SPE'] = np.mean(statsRec[traceID]['SPE'])
        statsRec[traceID]['median SPE'] = np.median(statsRec[traceID]['SPE'])

    for traceID in statsRec.keys():
        x.append(statsRec[traceID]['x'])
        y.append(statsRec[traceID]['y'])
        value.append(statsRec[traceID][variable])

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    size = []
    for val in value:
        size.append(100 * val / max(value))

    sc = ax.scatter(x, y, s=size, c=value)
    ax.text(0.5, 1.05,'Plot of all receivers',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    if not twoDim:
        ax.set_aspect('equal')
    cbar = ax.figure.colorbar(sc)
    cbar.set_label(variable)

    shot = survey.data.values()[0]
    for traceID in shot.getTraceIDlist():
        ax.annotate(' %s' % traceID, xy=(shot.getRecLoc(traceID)[0], shot.getRecLoc(traceID)[1]),
                    fontsize='x-small', color='k')
