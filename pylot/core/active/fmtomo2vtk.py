import numpy as np

def vgrids2VTK(inputfile = 'vgrids.in', outputfile = 'vgrids.vtk'):
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
        shotnumber = int(firstline.split()[1])
        nRayPoints = int(infile.readline().split()[0])
        if not shotnumber in rays.keys():
            rays[shotnumber] = {}
        rays[shotnumber][raynumber] = []
        for index in range(nRayPoints):
            if index%nthPoint is 0 or index == nRayPoints:
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

