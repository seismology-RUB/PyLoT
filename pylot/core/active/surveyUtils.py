import numpy as np

def generateSurvey(obsdir, shotlist):
    from obspy.core import read
    from pylot.core.active import seismicshot

    shot_dict = {}
    for shotnumber in shotlist:       # loop over data files
        # generate filenames and read manual picks to a list
        obsfile = obsdir + str(shotnumber) + '_pickle.dat'
        #obsfile = obsdir + str(shotnumber) + '.dat'
        
        if not obsfile in shot_dict.keys():
            shot_dict[shotnumber] = []
        shot_dict[shotnumber] = seismicshot.SeismicShot(obsfile)
        shot_dict[shotnumber].setParameters('shotnumber', shotnumber)
        
    return shot_dict

def setParametersForShots(cutwindow, tmovwind, tsignal, tgap, receiverfile, sourcefile, shot_dict):
    for shot in shot_dict.values():
        shot.setCut(cutwindow)
        shot.setTmovwind(tmovwind)
        shot.setTsignal(tsignal)
        shot.setTgap(tgap)
        shot.setRecfile(receiverfile)
        shot.setSourcefile(sourcefile)
        shot.setOrder(order = 4)

def removeEmptyTraces(shot_dict):
    filename = 'removeEmptyTraces.out'
    filename2 = 'updateTraces.out'
    outfile = open(filename, 'w')
    outfile2 = open(filename2, 'w')
    for shot in shot_dict.values():
        del_traceIDs = shot.updateTraceList()
        removed = shot.removeEmptyTraces()
        if removed is not None:
            outfile.writelines('shot: %s, removed empty traces: %s\n' %(shot.getShotnumber(), removed))
            outfile2.writelines('shot: %s, removed traceID(s) %s because they were not found in the corresponding stream\n' %(shot.getShotnumber(), del_traceIDs))
    print '\nremoveEmptyTraces, updateTraces: Finished! See %s and %s for more information of removed traces.\n' %(filename, filename2)
    outfile.close()
    outfile2.close()


def readParameters(parfile, parameter):
    from ConfigParser import ConfigParser
    parameterConfig = ConfigParser()
    parameterConfig.read('parfile')
    
    value = parameterConfig.get('vars', parameter).split('#')[0]
    value = value.replace(" ", "")

    return value

def setArtificialPick(shot_dict, traceID, pick):
    for shot in shot_dict.values():
        shot.setPick(traceID, pick)
        shot.setPickwindow(traceID, shot.getCut())

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

#def linearInterp(dist_med, dist_start

def exportFMTOMO(shot_dict, directory = 'FMTOMO_export', sourcefile = 'input_sf.in', ttFileExtension = '.tt'):
    count = 0
    fmtomo_factor = 1000 # transforming [m/s] -> [km/s]
    LatAll = []; LonAll = []; DepthAll = []
    srcfile = open(directory + '/' + sourcefile, 'w')
    srcfile.writelines('%10s\n' %len(shot_dict)) # number of sources
    for shotnumber in getShotlist(shot_dict):
        shot = getShotForShotnumber(shot_dict, shotnumber)
        ttfilename = str(shotnumber) + ttFileExtension
        (x, y, z) = shot.getSrcLoc() # getSrcLoc returns (x, y, z)
        srcfile.writelines('%10s %10s %10s\n' %(getAngle(y), getAngle(x), (-1)*z)) # lat, lon, depth
        LatAll.append(getAngle(y)); LonAll.append(getAngle(x)); DepthAll.append((-1)*z)
        srcfile.writelines('%10s\n' %1) # 
        srcfile.writelines('%10s %10s %10s\n' %(1, 1, ttfilename))
        ttfile = open(directory + '/' + ttfilename, 'w')
        traceIDlist = shot.getTraceIDlist()
        traceIDlist.sort()
        ttfile.writelines(str(countPickedTraces(shot)) + '\n')
        for traceID in traceIDlist:
            if shot.getPick(traceID) is not None:
                pick = shot.getPick(traceID) * fmtomo_factor
                delta = shot.getPickError(traceID) * fmtomo_factor
                (x, y, z) = shot.getRecLoc(traceID)
                ttfile.writelines('%20s %20s %20s %10s %10s\n' %(getAngle(y), getAngle(x), (-1)*z, pick, delta))
                LatAll.append(getAngle(y)); LonAll.append(getAngle(x)); DepthAll.append((-1)*z)
                count += 1
        ttfile.close()
    srcfile.close()       
    print 'Wrote output for %s traces' %count
    print 'WARNING: output generated for FMTOMO-obsdata. Obsdata seems to take Lat, Lon, Depth and creates output for FMTOMO as Depth, Lat, Lon'
    print 'Dimensions of the seismic Array, transformed for FMTOMO, are Depth(%s, %s), Lat(%s, %s), Lon(%s, %s)'%(
        min(DepthAll), max(DepthAll), min(LatAll), max(LatAll), min(LonAll), max(LonAll))

def getShotlist(shot_dict):
    shotlist = []
    for shot in shot_dict.values():
        shotlist.append(shot.getShotnumber())
    shotlist.sort()
    return shotlist

def getShotForShotnumber(shot_dict, shotnumber):
    for shot in shot_dict.values():
        if shot.getShotnumber() == shotnumber:
            return shot

def getAngle(distance):
    '''
    Function returns the angle on a Sphere of the radius R = 6371 [km] for a distance [km].
    '''
    PI = np.pi
    R = 6371.
    angle = distance * 180 / (PI * R)
    return angle

def countPickedTraces(shot):
    numtraces = 0
    for traceID in shot.getTraceIDlist():
        if shot.getPick(traceID) is not None:
            numtraces += 1
    print "countPickedTraces: Found %s picked traces in shot number %s"  %(numtraces, shot.getShotnumber())
    return numtraces

def countAllPickedTraces(shot_dict):
    traces = 0
    for shot in shot_dict.values():
        traces += countPickedTraces(shot)
    return traces

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
        
def vgrids2VTK(inputfile = 'vgrids.in', outputfile = 'vgrids.vtk'):
    def getDistance(angle):
        PI = np.pi
        R = 6371.
        distance = angle / 180 * (PI * R)
        return distance

    def readNumberOfPoints(filename):
        fin = open(filename, 'r')
        vglines = fin.readlines()

        nR = int(vglines[1].split()[0])
        nTheta = int(vglines[1].split()[1])
        nPhi = int(vglines[1].split()[2])

        fin.close()
        return nR, nTheta, nPhi

    def readDelta(filename):
        fin = open(filename, 'r')
        vglines = fin.readlines()

        dR = float(vglines[2].split()[0])
        dTheta = float(vglines[2].split()[1])
        dPhi = float(vglines[2].split()[2])

        fin.close()
        return dR, dTheta, dPhi

    def readStartpoints(filename):
        fin = open(filename, 'r')
        vglines = fin.readlines()

        sR = float(vglines[3].split()[0])
        sTheta = float(vglines[3].split()[1])
        sPhi = float(vglines[3].split()[2])

        fin.close()
        return sR, sTheta, sPhi

    def readVelocity(filename):
        vel = []; count = 0
        fin = open(filename, 'r')
        vglines = fin.readlines()

        for line in vglines:
            count += 1
            if count > 4:
                vel.append(float(line.split()[0]))

        print("Read %d points out of file: %s" %(count - 4, filename))
        return vel

    R = 6371 # earth radius
    outfile = open(outputfile, 'w')

    # Theta, Phi in radians, R in km
    nR, nTheta, nPhi = readNumberOfPoints(inputfile)
    dR, dTheta, dPhi = readDelta(inputfile)
    sR, sTheta, sPhi = readStartpoints(inputfile)
    vel = readVelocity(inputfile)

    nX = nPhi; nY = nTheta; nZ = nR

    sZ = sR - R
    sX = getDistance(np.rad2deg(sPhi))
    sY = getDistance(np.rad2deg(sTheta))

    dX = getDistance(np.rad2deg(dPhi))
    dY = getDistance(np.rad2deg(dTheta))
    
    xGrid = np.linspace(sX, sX + (dX * nX), nX)
    yGrid = np.linspace(sZ, sZ + (nY * dY), nY)
    zGrid = np.linspace(sZ, sZ + (nR * dR), nR)
    nPoints = len(xGrid) * len(yGrid) * len(zGrid)

    # write header
    print("Writing header for VTK file...")
    outfile.writelines('# vtk DataFile Version 3.1\n')
    outfile.writelines('Velocity on FMTOMO vgrids.in points\n')
    outfile.writelines('ASCII\n')
    outfile.writelines('DATASET POLYDATA\n')
    outfile.writelines('POINTS %15d float\n' %(nPoints))

    # write coordinates
    print("Writing coordinates to VTK file...")
    for z in zGrid:
        for y in yGrid:
            for x in xGrid:
                outfile.writelines('%10f %10f %10f \n' %(x, y, z))


    outfile.writelines('VERTICES %15d %15d\n' %(nPoints, 2 * nPoints))

    # write indices
    print("Writing indices to VTK file...")
    for index in range(nPoints):
        outfile.writelines('%10d %10d\n' %(1, index))

    outfile.writelines('POINT_DATA %15d\n' %(nPoints))
    outfile.writelines('SCALARS velocity float %d\n' %(1))
    outfile.writelines('LOOKUP_TABLE default\n')

    # write velocity
    print("Writing velocity values to VTK file...")
    for velocity in vel:
        outfile.writelines('%10f\n' %velocity)

    outfile.close()
    print("Wrote velocity grid for %d points to file: %s" %(nPoints, outputfile))
    return

