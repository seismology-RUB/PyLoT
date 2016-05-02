# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import datetime
import numpy as np

class Tomo3d(object):
    def __init__(self, nproc):
        self.defParas()
        self.nproc = nproc
        self.sources = self.readSrcFile()
        self.traces = self.readTraces()

    def defParas(self):
        self.fmm = 'fm3d'
        self.cvg = 'vgrids.in'
        self.cig = 'interfaces.in'
        self.csl = 'sources.in'
        self.pg = 'propgrid.in'
        self.rec = 'receivers.in'
        #self.ot = 'otimes.dat'
        self.frech = 'frechet.in'
        self.mode = 'mode_set.in'
        self.cwd = subprocess.check_output(['pwd'])[0:-1] + '/'
        self.folder = 'test_'

    def runFmm(self, directory, logfile, processes):
        os.chdir(directory)
        processes.append(subprocess.Popen('fm3d', stdout = None))
#        os.system('%s > %s &'%(self.fmm, logfile))
        os.chdir(self.cwd)
        return processes

    def makeDIR(self, directory):
        err = os.system('mkdir %s'%directory)
        if err is 256:
            response = raw_input('Warning: Directory already existing. Continue (y/n)?\n')
            if response == 'y':
                print('Overwriting existing files.')
            else:
                sys.exit('Aborted')
        
    def readNsrc(self):
        srcfile = open(self.csl, 'r')
        nsrc = int(srcfile.readline())
        srcfile.close()
        return nsrc

    def readNtraces(self):
        recfile = open(self.rec, 'r')
        nrec = int(recfile.readline())
        recfile.close()
        return nrec

    def calcSrcPerKernel(self):
        nsrc = self.readNsrc()
        return nsrc/self.nproc, nsrc%self.nproc

    def srcIDs4Kernel(self, procID):
        proc = procID - 1
        nsrc = self.readNsrc()
        srcPK, remain = self.calcSrcPerKernel()
        if procID > self.nproc:
            sys.exit('STOP: Kernel ID exceeds available number.')
        if proc < remain:
            start = (srcPK + 1) * (proc) + 1
            return range(start, start + srcPK + 1)
        elif proc == remain:
            start = (srcPK + 1) * (proc) + 1
            return range(start, start + srcPK)
        elif proc > remain:
            start = (srcPK + 1) * remain + srcPK * (proc - remain) + 1
            return range(start, start + srcPK)

    def startTomo(self):
        starttime = datetime.datetime.now()
        processes = []
        for procID in range(1, self.nproc + 1):
            directory = self.getProcDir(procID)
            log_out = self.cwd + 'fm3dlog_' + str(procID) + '.out'
            
            self.makeDIR(directory) # Problem bei Iteration
            self.writeSrcFile(procID, directory)
            self.writeTracesFile(procID, directory)
            os.system('cp %s %s %s %s %s %s %s' 
                      %(self.cvg, self.cig, self.frech,
                        self.fmm, self.mode, self.pg, directory))
            processes = self.runFmm(directory, log_out, processes)
            
        for p in processes:
            p.wait()

        self.mergeOutput()

        tdelta = datetime.datetime.now() - starttime
        print('Finished after %s'%tdelta)

    def readSrcFile(self):
        nsrc = self.readNsrc()
        srcfile = open(self.csl, 'r')

        sources = {}

        temp = srcfile.readline()
        for index in range(nsrc):
            teleflag = int(srcfile.readline())
            coords = srcfile.readline().split()
            numpaths = int(srcfile.readline())
            steps = int(srcfile.readline())
            interactions = srcfile.readline().split()
            veltype = int(srcfile.readline())
            if teleflag is not 0:
                sys.exit('Script not yet usable for teleseismic.')
            if numpaths is not 1:
                sys.exit('Script not yet usable for more than one path per source.')
            
            sources[index + 1] = {'teleflag': teleflag,
                                  'coords': coords,
                                  'numpaths': numpaths,
                                  'steps': steps,
                                  'interactions': interactions,
                                  'veltype': veltype
                                  }
            
        return sources

    def readTraces(self):
        recfile = open(self.rec, 'r')
        ntraces = self.readNtraces()
        
        traces = {}

        temp = recfile.readline()
        for index in range(ntraces):
            coords = recfile.readline().split()
            paths = int(recfile.readline())
            source = int(recfile.readline())
            path = int(recfile.readline())
            
            traces[index + 1] = { 'coords': coords,
                                  'paths': paths,
                                  'source': source,
                                  'path': path
                                  }
            
        return traces

    def writeSrcFile(self, procID, directory):
        srcfile = open('%s/sources.in'%directory, 'w')
        sourceIDs = self.srcIDs4Kernel(procID)

        srcfile.writelines('%s\n'%len(sourceIDs))
        for sourceID in sourceIDs:
            source = self.sources[sourceID]
            coords = source['coords']
            interactions = source['interactions']
            srcfile.writelines('%s\n'%source['teleflag'])
            srcfile.writelines('%s %s %s\n'%(float(coords[0]), float(coords[1]), float(coords[2])))
            srcfile.writelines('%s\n'%source['numpaths'])
            srcfile.writelines('%s\n'%source['steps'])
            srcfile.writelines('%s %s\n'%(int(interactions[0]), int(interactions[1])))
            srcfile.writelines('%s\n'%source['veltype'])

    def writeTracesFile(self, procID, directory):
        recfile = open('%s/receivers.in'%directory, 'w')
        sourceIDs = self.srcIDs4Kernel(procID)
        traceIDs = self.getTraceIDs4Sources(sourceIDs)

        recfile.writelines('%s\n'%len(traceIDs))
        for traceID in traceIDs:
            trace = self.traces[traceID]
            coords = trace['coords']
            source = int(trace['source']) - sourceIDs[0] + 1
            recfile.writelines('%s %s %s\n'%(float(coords[0]), float(coords[1]), float(coords[2])))
            recfile.writelines('%s\n'%trace['paths'])
            recfile.writelines('%s\n'%source)
            recfile.writelines('%s\n'%trace['path'])

    def getTraceIDs4Sources(self, sourceIDs):
        traceIDs = []
        for traceID in self.traces.keys():
            if self.traces[traceID]['source'] in sourceIDs:
                traceIDs.append(traceID)
        return traceIDs

    def getTraceIDs4Source(self, sourceID):
        traceIDs = []
        for traceID in self.traces.keys():
            if self.traces[traceID]['source'] == sourceID:
                traceIDs.append(traceID)
        return traceIDs

    def readArrivals(self, procID):
        directory = self.getProcDir(procID)
        arrfile = open(directory + '/arrivals.dat', 'r')
        sourceIDs = self.srcIDs4Kernel(procID)
        
        arrivals = []
        for sourceID in sourceIDs:
            traceIDs = self.getTraceIDs4Source(sourceID)
            for traceID in traceIDs:
                line = arrfile.readline().split()
                if line != []:
                    # recID and srcID for the individual processor will not be needed
                    recID_proc, srcID_proc, ray, normal, arrtime, diff, head = line
                    arrivals.append([traceID, sourceID, ray, normal, arrtime, diff, head])
                
        return arrivals

    def readFrechet(self, procID):
        directory = self.getProcDir(procID)
        frechfile = open(directory + '/frechet.dat', 'r')
        sourceIDs = self.srcIDs4Kernel(procID)
        
        frechet = {}
        for sourceID in sourceIDs:
            traceIDs = self.getTraceIDs4Source(sourceID)
            for traceID in traceIDs:
                line = frechfile.readline().split()
                if line != []:
                    # recID and srcID for the individual processor will not be needed
                    PDEV = []
                    recID_proc, srcID_proc, ray, normal, NPDEV = line
                    for i in range(int(NPDEV)):
                        PDEV.append(frechfile.readline())

                    frechet[traceID] = {'sourceID': sourceID,
                                        'raypath': ray,
                                        'normal': normal,
                                        'NPDEV': NPDEV,
                                        'PDEV': PDEV
                                        }
        return frechet

    def mergeArrivals(self):
        arrivalsOut = open(self.cwd + '/arrivals.dat', 'w')
        print('Merging arrivals.dat...')
        for procID in range(1, self.nproc + 1):
            arrivals = self.readArrivals(procID)
            for line in arrivals:
                arrivalsOut.writelines('%6s %6s %6s %6s %15s %5s %5s\n'%tuple(line))

    def mergeFrechet(self):
        print('Merging frechet.dat...')
        frechetOut = open(self.cwd + '/frechet.dat', 'w')
        for procID in range(1, self.nproc + 1):
            frechet = self.readFrechet(procID)
            traceIDs = frechet.keys()
            traceIDs.sort()
            for traceID in traceIDs:
                frech = frechet[traceID]
                frechetOut.writelines('%6s %6s %6s %6s %6s\n'%
                                      (traceID,
                                       frech['sourceID'],
                                       frech['raypath'],
                                       frech['normal'],
                                       frech['NPDEV']))
                for pdev in frech['PDEV']:
                    frechetOut.writelines(pdev)

    def mergeRays(self):
        print('Merging rays.dat...')
        filenames = []
        for procID in range(1, self.nproc + 1):
            directory = self.getProcDir(procID)
            filenames.append(directory + '/rays.dat')

        with open(self.cwd + 'rays.dat', 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

    def getProcDir(self, procID):
        return self.cwd + self.folder + str(procID)

    def mergeOutput(self):
        self.mergeArrivals()
        self.mergeFrechet()
        self.mergeRays()
        
        


def vgrids2VTK(inputfile='vgrids.in', outputfile='vgrids.vtk', absOrRel='abs', inputfileref='vgridsref.in'):
    '''
    Generate a vtk-file readable by e.g. paraview from FMTOMO output vgrids.in
    '''
    R = 6371.  # earth radius
    outfile = open(outputfile, 'w')

    number, delta, start, vel = _readVgrid(inputfile)

    nR, nTheta, nPhi = number
    dR, dTheta, dPhi = delta
    sR, sTheta, sPhi = start

    thetaGrid, phiGrid, rGrid = _generateGrids(number, delta, start)

    nPoints = nR * nTheta * nPhi

    nX = nPhi;
    nY = nTheta;
    nZ = nR

    sZ = sR - R
    sX = _getDistance(sPhi)
    sY = _getDistance(sTheta)

    dX = _getDistance(dPhi)
    dY = _getDistance(dTheta)
    dZ = dR

    # write header
    print("Writing header for VTK file...")
    outfile.writelines('# vtk DataFile Version 3.1\n')
    outfile.writelines('Velocity on FMTOMO vgrids.in points\n')
    outfile.writelines('ASCII\n')
    outfile.writelines('DATASET STRUCTURED_POINTS\n')

    outfile.writelines('DIMENSIONS %d %d %d\n' % (nX, nY, nZ))
    outfile.writelines('ORIGIN %f %f %f\n' % (sX, sY, sZ))
    outfile.writelines('SPACING %f %f %f\n' % (dX, dY, dZ))

    outfile.writelines('POINT_DATA %15d\n' % (nPoints))
    if absOrRel == 'abs':
<<<<<<< HEAD
        outfile.writelines('SCALARS velocity float %d\n' %(1))
    if absOrRel == 'relDepth':
        outfile.writelines('SCALARS velocity2depthMean float %d\n' %(1))
=======
        outfile.writelines('SCALARS velocity float %d\n' % (1))
>>>>>>> 37f9292c39246b327d3630995ca2521725c6cdd7
    elif absOrRel == 'rel':
        outfile.writelines('SCALARS velChangePercent float %d\n' % (1))
    outfile.writelines('LOOKUP_TABLE default\n')

    pointsPerR = nTheta * nPhi

    # write velocity
    if absOrRel == 'abs':
        print("Writing velocity values to VTK file...")
        for velocity in vel:
<<<<<<< HEAD
            outfile.writelines('%10f\n' %velocity)
    elif absOrRel == 'relDepth':
        print("Writing velocity values to VTK file relative to mean of each depth...")
        index = 0; count = 0
        veldepth = []
        for velocity in vel:
            count += 1
            veldepth.append(velocity)
            if count%pointsPerR == 0:
                velmean = np.mean(veldepth)
                #print velmean, count, count/pointsPerR
                for vel in veldepth:
                    outfile.writelines('%10f\n' %(vel - velmean))
                veldepth = []
=======
            outfile.writelines('%10f\n' % velocity)
>>>>>>> 37f9292c39246b327d3630995ca2521725c6cdd7
    elif absOrRel == 'rel':
        nref, dref, sref, velref = _readVgrid(inputfileref)
        nR_ref, nTheta_ref, nPhi_ref = nref
        if not len(velref) == len(vel):
            print('ERROR: Number of gridpoints mismatch for %s and %s' % (inputfile, inputfileref))
            return
        # velrel = [((vel - velref) / velref * 100) for vel, velref in zip(vel, velref)]
        velrel = []
        for velocities in zip(vel, velref):
            v, vref = velocities
            if not vref == 0:
                velrel.append((v - vref) / vref * 100)
            else:
                velrel.append(0)

        if not nR_ref == nR and nTheta_ref == nTheta and nPhi_ref == nPhi:
            print('ERROR: Dimension mismatch of grids %s and %s' % (inputfile, inputfileref))
            return
        print("Writing velocity values to VTK file...")
        for velocity in velrel:
            outfile.writelines('%10f\n' % velocity)
        print('Pertubations: min: %s %%, max: %s %%' % (min(velrel), max(velrel)))

    outfile.close()
    print("Wrote velocity grid for %d points to file: %s" % (nPoints, outputfile))
    return


def rays2VTK(fnin, fdirout='./vtk_files/', nthPoint=50):
    '''
    Writes VTK file(s) for FMTOMO rays from rays.dat

    :param: nthPoint, plot every nth point of the ray
    :type: integer
    '''
    infile = open(fnin, 'r')
    R = 6371
    rays = {}
    raynumber = 0
    nPoints = 0

    ### NOTE: rays.dat seems to be in km and radians
    while True:
        raynumber += 1
        firstline = infile.readline()
        if firstline == '': break  # break at EOF
        raynumber = int(firstline.split()[0])
        shotnumber = int(firstline.split()[1])
        rayValid = int(firstline.split()[4])  # is zero if the ray is invalid
        if rayValid == 0:
            print('Invalid ray number %d for shot number %d' % (raynumber, shotnumber))
            continue
        nRayPoints = int(infile.readline().split()[0])
        if not shotnumber in rays.keys():
            rays[shotnumber] = {}
        rays[shotnumber][raynumber] = []
        for index in range(nRayPoints):
            if index % nthPoint is 0 or index == (nRayPoints - 1):
                rad, lat, lon = infile.readline().split()
                rays[shotnumber][raynumber].append(
                    [_getDistance(np.rad2deg(float(lon))), _getDistance(np.rad2deg(float(lat))), float(rad) - R])
            else:
                dummy = infile.readline()

    infile.close()

    for shotnumber in rays.keys():
        fnameout = fdirout + 'rays%03d.vtk' % (shotnumber)
        outfile = open(fnameout, 'w')

        nPoints = 0
        for raynumber in rays[shotnumber]:
            for ray in rays[shotnumber][raynumber]:
                nPoints += 1

        # write header
        # print("Writing header for VTK file...")
        print("Writing shot %d to file %s" % (shotnumber, fnameout))
        outfile.writelines('# vtk DataFile Version 3.1\n')
        outfile.writelines('FMTOMO rays\n')
        outfile.writelines('ASCII\n')
        outfile.writelines('DATASET POLYDATA\n')
        outfile.writelines('POINTS %15d float\n' % (nPoints))

        # write coordinates
        # print("Writing coordinates to VTK file...")
        for raynumber in rays[shotnumber].keys():
            for raypoint in rays[shotnumber][raynumber]:
                outfile.writelines('%10f %10f %10f \n' % (raypoint[0], raypoint[1], raypoint[2]))

        outfile.writelines('LINES %15d %15d\n' % (len(rays[shotnumber]), len(rays[shotnumber]) + nPoints))

        # write indices
        # print("Writing indices to VTK file...")
        count = 0
        for raynumber in rays[shotnumber].keys():
            outfile.writelines('%d ' % (len(rays[shotnumber][raynumber])))
            for index in range(len(rays[shotnumber][raynumber])):
                outfile.writelines('%d ' % (count))
                count += 1
            outfile.writelines('\n')


def _readVgrid(filename):
    def readNumberOfPoints(filename):
        fin = open(filename, 'r')
        vglines = fin.readlines()

        nR = int(vglines[1].split()[0])
        nTheta = int(vglines[1].split()[1])
        nPhi = int(vglines[1].split()[2])

        print('readNumberOf Points: Awaiting %d grid points in %s'
              % (nR * nTheta * nPhi, filename))
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
        '''
        Reads in velocity from vgrids file and returns a list containing all values in the same order
        '''
        vel = [];
        count = 0
        fin = open(filename, 'r')
        vglines = fin.readlines()

        for line in vglines:
            count += 1
            if count > 4:
                vel.append(float(line.split()[0]))

        print("Read %d points out of file: %s" % (count - 4, filename))
        return vel

    # Theta, Phi in radians, R in km
    nR, nTheta, nPhi = readNumberOfPoints(filename)
    dR, dThetaRad, dPhiRad = readDelta(filename)
    sR, sThetaRad, sPhiRad = readStartpoints(filename)
    vel = readVelocity(filename)

    dTheta, dPhi = np.rad2deg((dThetaRad, dPhiRad))
    sTheta, sPhi = np.rad2deg((sThetaRad, sPhiRad))

    number = (nR, nTheta, nPhi)
    delta = (dR, dTheta, dPhi)
    start = (sR, sTheta, sPhi)
    return number, delta, start, vel


def _generateGrids(number, delta, start):
    nR, nTheta, nPhi = number
    dR, dTheta, dPhi = delta
    sR, sTheta, sPhi = start

    eR = sR + (nR - 1) * dR
    ePhi = sPhi + (nPhi - 1) * dPhi
    eTheta = sTheta + (nTheta - 1) * dTheta

    thetaGrid = np.linspace(sTheta, eTheta, num=nTheta)
    phiGrid = np.linspace(sPhi, ePhi, num=nPhi)
    rGrid = np.linspace(sR, eR, num=nR)

    return (thetaGrid, phiGrid, rGrid)


def addCheckerboard(spacing=10., pertubation=0.1, inputfile='vgrids.in',
                    outputfile='vgrids_cb.in', ampmethod='linear', rect=(None, None)):
    '''
    Add a checkerboard to an existing vgrids.in velocity model.

    :param: spacing, size of the tiles
    type: float

    :param: pertubation, pertubation (default: 0.1 = 10%)
    type: float
    '''

    def correctSpacing(spacing, delta, disttype=None):
        if spacing > delta:
            spacing_corr = round(spacing / delta) * delta
        elif spacing < delta:
            spacing_corr = delta
        print('The spacing of the checkerboard of %s (%s) was corrected to '
              'a value of %s to fit the grid spacing of %s.' % (spacing, disttype, spacing_corr, delta))
        return spacing_corr

    def linearAmp(InCell):
        decimal = InCell - np.floor(InCell)
        return (-abs(decimal - 0.5) + 0.5) * 2

    def rectAmp(InCell, rect):
        decimal = InCell - np.floor(InCell)
        r1, r2 = rect
        if r1 <= decimal <= r2:
            return 1
        else:
            return 0

    def ampFunc(InCell, method='linear', rect=None):
        if method == 'linear':
            return linearAmp(InCell)
        if method == 'rect' and rect is not None:
            return rectAmp(InCell, rect)
        else:
            print('ampFunc: Could not amplify cb pattern')

    decm = 0.3  # diagonal elements of the covariance matrix (grid3dg's default value is 0.3)
    outfile = open(outputfile, 'w')

    number, delta, start, vel = _readVgrid(inputfile)

    nR, nTheta, nPhi = number
    dR, dTheta, dPhi = delta
    sR, sTheta, sPhi = start

    thetaGrid, phiGrid, rGrid = _generateGrids(number, delta, start)

    nPoints = nR * nTheta * nPhi

    # write header for velocity grid file (in RADIANS)
    outfile.writelines('%10s %10s \n' % (1, 1))
    outfile.writelines('%10s %10s %10s\n' % (nR, nTheta, nPhi))
    outfile.writelines('%10s %10s %10s\n' % (dR, np.deg2rad(dTheta), np.deg2rad(dPhi)))
    outfile.writelines('%10s %10s %10s\n' % (sR, np.deg2rad(sTheta), np.deg2rad(sPhi)))

    spacR = correctSpacing(spacing, dR, '[meter], R')
    spacTheta = correctSpacing(_getAngle(spacing), dTheta, '[degree], Theta')
    spacPhi = correctSpacing(_getAngle(spacing), dPhi, '[degree], Phi')

    count = 0
    evenOdd = 1
    even = 0;
    odd = 0

    # In the following loop it is checked whether the positive distance from the border of the model
    # for a point on the grid divided by the spacing is even or odd and then pertubated.
    # The position is also shifted by half of the delta so that the position is directly on the point and
    # not on the border between two points.
    # "InCell" points e.g. rInCell are floats with their integer number corresponding to the cell number and
    # their decimal place (0 - 1) corresponding to the position inside the cell.
    # The amplification factor ampFactor comes from a linear relationship and ranges between 0 (cell border)
    # and 1 (cell middle)
    for radius in rGrid:
        rInCell = (radius - sR - dR / 2) / spacR
        ampR = ampFunc(rInCell, ampmethod, rect)
        if np.floor(rInCell) % 2:
            evenOddR = 1
        else:
            evenOddR = -1
        for theta in thetaGrid:
            thetaInCell = (theta - sTheta - dTheta / 2) / spacTheta
            ampTheta = ampFunc(thetaInCell, ampmethod, rect)
            if np.floor(thetaInCell) % 2:
                evenOddT = 1
            else:
                evenOddT = -1
            for phi in phiGrid:
                phiInCell = (phi - sPhi - dPhi / 2) / spacPhi
                ampPhi = ampFunc(phiInCell, ampmethod, rect)
                if np.floor(phiInCell) % 2:
                    evenOddP = 1
                else:
                    evenOddP = -1
                velocity = vel[count]
                ampFactor = (ampR + ampTheta + ampPhi) / 3
                evenOdd = evenOddR * evenOddT * evenOddP * ampFactor
                velocity += evenOdd * pertubation * velocity

                outfile.writelines('%10s %10s\n' % (velocity, decm))
                count += 1

                progress = float(count) / float(nPoints) * 100
                _update_progress(progress)

    print('Added checkerboard to the grid in file %s with a spacing of %s and a pertubation of %s %%. '
          'Outputfile: %s.' % (inputfile, spacing, pertubation * 100, outputfile))
    outfile.close()


def addBox(x=(None, None), y=(None, None), z=(None, None),
           boxvelocity=1.0, inputfile='vgrids.in',
           outputfile='vgrids_box.in'):
    '''
    Add a box with constant velocity to an existing vgrids.in velocity model.

    :param: x, borders of the box (xleft, xright)
    type: tuple

    :param: y, borders of the box (yleft, yright)
    type: tuple

    :param: z, borders of the box (bot, top)
    type: tuple

    :param: boxvelocity, default: 1.0 km/s
    type: float
    '''
    R = 6371.
    decm = 0.3  # diagonal elements of the covariance matrix (grid3dg's default value is 0.3)
    outfile = open(outputfile, 'w')

    theta1 = _getAngle(y[0])
    theta2 = _getAngle(y[1])
    phi1 = _getAngle(x[0])
    phi2 = _getAngle(x[1])
    r1 = R + z[0]
    r2 = R + z[1]

    print('Adding box to grid with theta = (%s, %s), phi = (%s, %s), '
          'r = (%s, %s), velocity = %s [km/s]'
          % (theta1, theta2, phi1, phi2, r1, r2, boxvelocity))

    number, delta, start, vel = _readVgrid(inputfile)

    nR, nTheta, nPhi = number
    dR, dTheta, dPhi = delta
    sR, sTheta, sPhi = start

    thetaGrid, phiGrid, rGrid = _generateGrids(number, delta, start)

    nPoints = nR * nTheta * nPhi

    # write header for velocity grid file (in RADIANS)
    outfile.writelines('%10s %10s \n' % (1, 1))
    outfile.writelines('%10s %10s %10s\n' % (nR, nTheta, nPhi))
    outfile.writelines('%10s %10s %10s\n' % (dR, np.deg2rad(dTheta), np.deg2rad(dPhi)))
    outfile.writelines('%10s %10s %10s\n' % (sR, np.deg2rad(sTheta), np.deg2rad(sPhi)))

    count = 0
    for radius in rGrid:
        if r1 <= radius <= r2:
            rFlag = 1
        else:
            rFlag = 0
        for theta in thetaGrid:
            if theta1 <= theta <= theta2:
                thetaFlag = 1
            else:
                thetaFlag = 0
            for phi in phiGrid:
                if phi1 <= phi <= phi2:
                    phiFlag = 1
                else:
                    phiFlag = 0
                velocity = vel[count]
                if rFlag * thetaFlag * phiFlag is not 0:
                    velocity = boxvelocity

                outfile.writelines('%10s %10s\n' % (velocity, decm))
                count += 1

                progress = float(count) / float(nPoints) * 100
                _update_progress(progress)

    print('Added box to the grid in file %s. '
          'Outputfile: %s.' % (inputfile, outputfile))
    outfile.close()


def _update_progress(progress):
    sys.stdout.write("%d%% done   \r" % (progress))
    sys.stdout.flush()


def _getAngle(distance):
    '''
    Function returns the angle on a Sphere of the radius R = 6371 [km] for a distance [km].
    '''
    PI = np.pi
    R = 6371.
    angle = distance * 180. / (PI * R)
    return angle


def _getDistance(angle):
    PI = np.pi
    R = 6371.
    distance = angle / 180 * (PI * R)
    return distance
