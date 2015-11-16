# -*- coding: utf-8 -*-
import sys
import numpy as np

def vgrids2VTK(inputfile = 'vgrids.in', outputfile = 'vgrids.vtk', absOrRel = 'abs', inputfileref = 'vgridsref.in'):
    '''
    Generate a vtk-file readable by e.g. paraview from FMTOMO output vgrids.in
    '''
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

        print('readNumberOf Points: Awaiting %d grid points in %s'
              %(nR*nTheta*nPhi, filename))
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
        vel = []; count = 0
        fin = open(filename, 'r')
        vglines = fin.readlines()

        for line in vglines:
            count += 1
            if count > 4:
                vel.append(float(line.split()[0]))

        print("Read %d points out of file: %s" %(count - 4, filename))
        return vel

    R = 6371. # earth radius
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

    nPoints = nX * nY * nZ

    dZ = dR

    # write header
    print("Writing header for VTK file...")
    outfile.writelines('# vtk DataFile Version 3.1\n')
    outfile.writelines('Velocity on FMTOMO vgrids.in points\n')
    outfile.writelines('ASCII\n')
    outfile.writelines('DATASET STRUCTURED_POINTS\n')

    outfile.writelines('DIMENSIONS %d %d %d\n' %(nX, nY, nZ))
    outfile.writelines('ORIGIN %f %f %f\n' %(sX, sY, sZ))
    outfile.writelines('SPACING %f %f %f\n' %(dX, dY, dZ))

    outfile.writelines('POINT_DATA %15d\n' %(nPoints))
    if absOrRel == 'abs':
        outfile.writelines('SCALARS velocity float %d\n' %(1))
    elif absOrRel == 'rel':
        outfile.writelines('SCALARS velChangePercent float %d\n' %(1))
    outfile.writelines('LOOKUP_TABLE default\n')

    # write velocity
    if absOrRel == 'abs':
        print("Writing velocity values to VTK file...")
        for velocity in vel:
            outfile.writelines('%10f\n' %velocity)
    elif absOrRel == 'rel':
        velref = readVelocity(inputfileref)
        if not len(velref) == len(vel):
            print('ERROR: Number of gridpoints mismatch for %s and %s'%(inputfile, inputfileref))
            return
        #velrel = [((vel - velref) / velref * 100) for vel, velref in zip(vel, velref)]
        velrel = []
        for velocities in zip(vel, velref):
            v, vref = velocities
            if not vref == 0:
                velrel.append((v - vref) / vref * 100)
            else:
                velrel.append(0)            

        nR_ref, nTheta_ref, nPhi_ref = readNumberOfPoints(inputfileref)
        if not nR_ref == nR and nTheta_ref == nTheta and nPhi_ref == nPhi:
            print('ERROR: Dimension mismatch of grids %s and %s'%(inputfile, inputfileref))
            return
        print("Writing velocity values to VTK file...")
        for velocity in velrel:
            outfile.writelines('%10f\n' %velocity)
        print('Pertubations: min: %s %%, max: %s %%'%(min(velrel), max(velrel)))

    outfile.close()
    print("Wrote velocity grid for %d points to file: %s" %(nPoints, outputfile))
    return

def rays2VTK(fnin, fdirout = './vtk_files/', nthPoint = 50):
    '''
    Writes VTK file(s) for FMTOMO rays from rays.dat

    :param: nthPoint, plot every nth point of the ray
    :type: integer
    '''
    def getDistance(angle):
        PI = np.pi
        R = 6371.
        distance = angle / 180 * (PI * R)
        return distance

    infile = open(fnin, 'r')
    R = 6371
    rays = {}
    raynumber = 0
    nPoints = 0

    ### NOTE: rays.dat seems to be in km and radians

    while True:
        raynumber += 1
        firstline = infile.readline()
        if firstline == '': break # break at EOF
        raynumber = int(firstline.split()[0])
        shotnumber = int(firstline.split()[1])
        rayValid = int(firstline.split()[4]) # is zero if the ray is invalid
        if rayValid == 0:
            print('Invalid ray number %d for shot number %d'%(raynumber, shotnumber))
            continue
        nRayPoints = int(infile.readline().split()[0])
        if not shotnumber in rays.keys():
            rays[shotnumber] = {}
        rays[shotnumber][raynumber] = []
        for index in range(nRayPoints):
            if index % nthPoint is 0 or index == (nRayPoints - 1):
                rad, lat, lon = infile.readline().split()
                rays[shotnumber][raynumber].append([getDistance(np.rad2deg(float(lon))), getDistance(np.rad2deg(float(lat))), float(rad) - R])
            else:
                dummy = infile.readline()

    infile.close()

    for shotnumber in rays.keys():
        fnameout = fdirout + 'rays%03d.vtk'%(shotnumber)
        outfile = open(fnameout, 'w')

        nPoints = 0
        for raynumber in rays[shotnumber]:
            for ray in rays[shotnumber][raynumber]:
                nPoints += 1

        # write header
        #print("Writing header for VTK file...")
        print("Writing shot %d to file %s" %(shotnumber, fnameout))
        outfile.writelines('# vtk DataFile Version 3.1\n')
        outfile.writelines('FMTOMO rays\n')
        outfile.writelines('ASCII\n')
        outfile.writelines('DATASET POLYDATA\n')
        outfile.writelines('POINTS %15d float\n' %(nPoints))

        # write coordinates
        #print("Writing coordinates to VTK file...")

        for raynumber in rays[shotnumber].keys():
            for raypoint in rays[shotnumber][raynumber]:
                outfile.writelines('%10f %10f %10f \n' %(raypoint[0], raypoint[1], raypoint[2]))

        outfile.writelines('LINES %15d %15d\n' %(len(rays[shotnumber]), len(rays[shotnumber]) + nPoints))

        # write indices
        #print("Writing indices to VTK file...")

        count = 0
        for raynumber in rays[shotnumber].keys():
            outfile.writelines('%d ' %(len(rays[shotnumber][raynumber])))
            for index in range(len(rays[shotnumber][raynumber])):
                outfile.writelines('%d ' %(count))
                count += 1
            outfile.writelines('\n')

    # outfile.writelines('POINT_DATA %15d\n' %(nPoints))
    # outfile.writelines('SCALARS rays float %d\n' %(1))
    # outfile.writelines('LOOKUP_TABLE default\n')

    # # write velocity
    # print("Writing velocity values to VTK file...")
    # for velocity in vel:
    #     outfile.writelines('%10f\n' %velocity)

    # outfile.close()
    # print("Wrote velocity grid for %d points to file: %s" %(nPoints, outputfile))

def addCheckerboard(spacing = 10., pertubation = 0.1, inputfile = 'vgrids.in',
                    outputfile = 'vgrids_cb.in', ampmethod = 'linear', rect = (None, None)):
    '''
    Add a checkerboard to an existing vgrids.in velocity model.

    :param: spacing, size of the tiles
    type: float

    :param: pertubation, pertubation (default: 0.1 = 10%)
    type: float
    '''
    def readNumberOfPoints(filename):
        fin = open(filename, 'r')
        vglines = fin.readlines()

        nR = int(vglines[1].split()[0])
        nTheta = int(vglines[1].split()[1])
        nPhi = int(vglines[1].split()[2])

        print('readNumberOf Points: Awaiting %d grid points in %s'
              %(nR*nTheta*nPhi, filename))
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
        vel = []; count = 0
        fin = open(filename, 'r')
        vglines = fin.readlines()

        for line in vglines:
            count += 1
            if count > 4:
                vel.append(float(line.split()[0]))

        print("Read %d points out of file: %s" %(count - 4, filename))
        return vel

    def correctSpacing(spacing, delta, disttype = None):
        if spacing > delta:
            spacing_corr = round(spacing / delta) * delta
        elif spacing < delta:
            spacing_corr = delta
        print('The spacing of the checkerboard of %s (%s) was corrected to '
              'a value of %s to fit the grid spacing of %s.' %(spacing, disttype, spacing_corr, delta))
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

    def ampFunc(InCell, method = 'linear', rect = None):
        if method == 'linear':
            return linearAmp(InCell)
        if method == 'rect' and rect is not None:
            return rectAmp(InCell, rect)
        else:
            print('ampFunc: Could not amplify cb pattern')


    R = 6371. # earth radius
    decm = 0.3 # diagonal elements of the covariance matrix (grid3dg's default value is 0.3)
    outfile = open(outputfile, 'w')

    # Theta, Phi in radians, R in km
    nR, nTheta, nPhi = readNumberOfPoints(inputfile)
    dR, dThetaRad, dPhiRad = readDelta(inputfile)
    sR, sThetaRad, sPhiRad = readStartpoints(inputfile)
    vel = readVelocity(inputfile)

    dTheta, dPhi = np.rad2deg((dThetaRad, dPhiRad))
    sTheta, sPhi = np.rad2deg((dThetaRad, dPhiRad))

    eR = sR + (nR - 1) * dR
    ePhi = sPhi + (nPhi - 1) * dPhi
    eTheta = sTheta + (nTheta - 1) * dTheta

    nPoints = nR * nTheta * nPhi

    thetaGrid = np.linspace(sTheta, eTheta, num = nTheta)
    phiGrid = np.linspace(sPhi, ePhi, num = nPhi)
    rGrid = np.linspace(sR, eR, num = nR)

    # write header for velocity grid file (in RADIANS)
    outfile.writelines('%10s %10s \n' %(1, 1))
    outfile.writelines('%10s %10s %10s\n' %(nR, nTheta, nPhi))
    outfile.writelines('%10s %10s %10s\n' %(dR, dThetaRad, dPhiRad))
    outfile.writelines('%10s %10s %10s\n' %(sR, sThetaRad, sPhiRad))

    spacR = correctSpacing(spacing, dR, '[meter], R')
    spacTheta = correctSpacing(_getAngle(spacing), dTheta, '[degree], Theta')
    spacPhi = correctSpacing(_getAngle(spacing), dPhi, '[degree], Phi')

    count = 0
    evenOdd = 1
    even = 0; odd = 0

    # In the following loop it is checked whether the positive distance from the border of the model
    # for a point on the grid divided by the spacing is even or odd and then pertubated.
    # The position is also shifted by half of the delta so that the position is directly on the point and
    # not on the border between two points.
    # "InCell" points e.g. rInCell are floats with their integer number corresponding to the cell number and
    # their decimal place (0 - 1) corresponding to the position inside the cell.
    # The amplification factor ampFactor comes from a linear relationship and ranges between 0 (cell border)
    # and 1 (cell middle)
    for radius in rGrid:
        rInCell = (radius - sR - dR/2) / spacR 
        ampR = ampFunc(rInCell, ampmethod, rect)
        if np.floor(rInCell) % 2:
            evenOddR = 1
        else:
            evenOddR = -1
        for theta in thetaGrid:
            thetaInCell = (theta - sTheta - dTheta/2) / spacTheta
            ampTheta = ampFunc(thetaInCell, ampmethod, rect)
            if np.floor(thetaInCell) % 2:
                evenOddT = 1
            else:
                evenOddT = -1
            for phi in phiGrid:
                phiInCell = (phi - sPhi - dPhi/2) / spacPhi
                ampPhi = ampFunc(phiInCell, ampmethod, rect)
                if np.floor(phiInCell) % 2:
                    evenOddP = 1
                else:
                    evenOddP = -1
                velocity = vel[count]
                ampFactor = (ampR + ampTheta + ampPhi) / 3
                evenOdd = evenOddR * evenOddT * evenOddP * ampFactor
                velocity += evenOdd * pertubation * velocity

                outfile.writelines('%10s %10s\n'%(velocity, decm))
                count += 1

                progress = float(count) / float(nPoints) * 100
                _update_progress(progress)

    print('Added checkerboard to the grid in file %s with a spacing of %s and a pertubation of %s %%. '
          'Outputfile: %s.'%(inputfile, spacing, pertubation, outputfile))
    outfile.close()

def _update_progress(progress):
    sys.stdout.write("%d%% done   \r" % (progress) )
    sys.stdout.flush()

def _getAngle(distance):
    '''
    Function returns the angle on a Sphere of the radius R = 6371 [km] for a distance [km].
    '''
    PI = np.pi
    R = 6371.
    angle = distance * 180. / (PI * R)
    return angle

