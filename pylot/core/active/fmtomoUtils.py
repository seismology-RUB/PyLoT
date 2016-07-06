# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import datetime
import numpy as np

class Tomo3d(object):
    def __init__(self, fmtomodir, simuldir = 'fmtomo_simulation', citer = 0, overwrite = False):
        '''
        Class build from FMTOMO script tomo3d. Can be used to run several instances of FMM code in parallel.

        :param: citer, current iteration (default = 0: start new model)
        :type: integer
        '''
        self.simuldir = simuldir
        self.setCWD()
        self.buildFmtomodir(fmtomodir)
        self.buildObsdata()
        self.defParas()
        self.copyRef()
        self.citer = citer       # current iteration
        self.sources = self.readSrcFile()
        self.traces = self.readTraces()
        self.directories = []
        self.overwrite = overwrite
        
    def defParas(self):
        self.defFMMParas()
        self.defInvParas()

    def buildFmtomodir(self, directory):
        tomo_files = ['fm3d',
        'frechgen',
        'frechgen.in',
        'invert3d',
        'invert3d.in',
        'mode_set.in',
        'obsdata',
        'obsdata.in',
        'residuals',
        'residuals.in',
        'tomo3d',
        'tomo3d.in']

        for name in tomo_files:
            filename = os.path.join(directory, name)
            linkname = os.path.join(self.cwd, name)
            os.system('ln -s %s %s'%(filename, linkname))

    def buildObsdata(self):
        os.system('obsdata')
        os.system('mv sources.in sourcesref.in')

    def defFMMParas(self):
        '''
        Initiates parameters for the forward calculation.
        '''
        # Name of fast marching program
        self.fmm = '{0}/fm3d'.format(self.cwd)
        # Name of program calculating Frechet derivatives
        self.frechgen = '{0}/frechgen'.format(self.cwd)
        # Name of current velocity/inversion grid
        self.cvg = 'vgrids.in'
        # Name of current interfaces grid
        self.cig = 'interfaces.in'
        # Name of file containing current source locations
        self.csl = 'sources.in'
        # Name of file containing propagation grid
        self.pg = 'propgrid.in'
        # Name of file containing receiver coordinates
        self.rec = 'receivers.in'
        self.frech = 'frechet.in'
        self.frechout = 'frechet.dat'
        # Name of file containing measured data
        self.ot = 'otimes.dat'
        # Name of file containing output velocity information
        self.ttim = 'arrivals.dat'
        self.mode = 'mode_set.in'
        # Name of temporary folders created for each process
        self.folder = '.proc_'

    def defInvParas(self):
        '''
        Initiates inversion parameters for FMTOMO.
        '''
        # Name of program for performing inversion
        self.inv = '{0}/invert3d'.format(self.cwd)
        # Name of file containing current model traveltimes
        self.mtrav = 'mtimes.dat'
        # Name of file containing reference model traveltimes
        self.rtrav = 'rtimes.dat'
        # Name of file containing initial velocity grid
        self.ivg = 'vgridsref.in'
        # Name of file containing initial interface grid
        self.iig = 'interfacesref.in'
        # Name of file containing initial source locations
        self.isl = 'sourcesref.in'
        # Name of program for calculating traveltime residuals
        self.resid = '{0}/residuals'.format(self.cwd)
        # Name of output file for calculating traveltime residuals
        self.resout = 'residuals.dat'

    def copyRef(self):
        '''
        Copies reference grids to used grids (e.g. sourcesref.in to sources.in)
        '''
        os.system('cp %s %s'%(self.ivg, self.cvg))
        os.system('cp %s %s'%(self.iig, self.cig))
        os.system('cp %s %s'%(self.isl, self.csl))

    def setCWD(self, directory = None):
        '''
        Set working directory containing all necessary files.

        Default: pwd
        '''
        if directory == None:
            directory = os.path.join(os.getcwd(), self.simuldir)

        os.chdir(directory)
        self.cwd = directory
        print('Working directory is: %s'%self.cwd)

    def runFrech(self):
        os.system(self.frechgen)

    def runTOMO3D(self, nproc, iterations):
        '''
        Starts up the FMTOMO code for the set number of iterations on nproc parallel processes.
        
        :param: nproc, number of parallel processes
        :type: integer

        :param: iterations, number of iterations
        :type: integer
        '''
        self.nproc = nproc
        self.iter = iterations   # number of iterations

        starttime = datetime.datetime.now()
        print('Starting TOMO3D on %s parallel processes for %s iteration(s).'
              %(self.nproc, self.iter))
        if self.citer == 0:
            self.makeInvIterDir()
            self.startForward(self.cInvIterDir)
            self.raiseIter()
        
        while self.citer <= self.iter:
            self.makeInvIterDir()
            self.startInversion()
            self.saveVgrid()
            self.startForward(self.cInvIterDir)
            self.raiseIter()

        if self.citer > self.iter:
            self.removeDirectories()
            self.unlink(os.path.join(self.cwd, self.frechout))
            self.unlink(os.path.join(self.cwd, self.ttim))

        tdelta = datetime.datetime.now() - starttime
        print('runTOMO3D: Finished %s iterations after %s.'%(self.iter, tdelta))

    def runFmm(self, directory, logfile, processes):
        '''
        Calls an instance of the FMM code in the process directory.
        Requires a list of all active processes and returns an updated list.
        '''
        os.chdir(directory)
        fout = open(logfile, 'w')
        processes.append(subprocess.Popen(self.fmm, stdout = fout))
        fout.close()
        os.chdir(self.cwd)
        return processes

    def startForward(self, logdir):
        '''
        Runs an instance of the FMM code in the process directory.
        '''
        self._printLine()
        print('Starting forward simulation for iteration %s.'%(self.citer))

        if self.citer == 0:
            self.copyRef()
            self.runFrech()
            self.makeDirectories()

        starttime = datetime.datetime.now()
        processes = []

        for procID in range(1, self.nproc + 1):
            directory = self.getProcDir(procID)
            logfn = 'fm3dlog_' + str(procID) + '.out'
            log_out = os.path.join(logdir, logfn)
            
            self.writeSrcFile(procID)
            self.writeTracesFile(procID)
            os.system('cp {cvg} {cig} {mode} {pg} {frechin} {dest}'
                      .format(cvg=self.cvg, cig=self.cig, frechin=self.frech,
                              mode=self.mode, pg=self.pg, dest=directory))
            processes = self.runFmm(directory, log_out, processes)
            
        for p in processes:
            p.wait()

        self.mergeOutput(self.cInvIterDir)
        self.clearDirectories()
        self.copyArrivals()
        if self.citer == 0:
            self.copyArrivals(self.rtrav)

        self.calcRes()
        tdelta = datetime.datetime.now() - starttime
        print('Finished Forward calculation after %s'%tdelta)

    def startInversion(self):
        '''
        Simply calls the inversion program.
        '''
        print('Calling %s...'%self.inv)
        os.system(self.inv)

    def calcRes(self):
        '''
        Calls residual calculation program.
        '''
        resout = os.path.join(self.cwd, self.resout)
        if self.citer == 0:
            os.system('%s > %s'%(self.resid, resout))
        else:
            os.system('%s >> %s'%(self.resid, resout))

        with open(resout, 'r') as infile:
            residuals = infile.readlines()
        RMS, var, chi2 = residuals[-1].split()
        print('Residuals: RMS = %s, var = %s, Chi^2 = %s.'%(RMS, var, chi2))

    def raiseIter(self):
        self.citer +=1
        self._printLine()
        invfile = open(self.cwd + '/inviter.in', 'w')
        invfile.write('%s'%self.citer)
        invfile.close()

    def makeDir(self, directory):
        err = os.system('mkdir %s'%directory)
        if err is 0:
            self.directories.append(directory)
            return
        if err is 256:
            if self.overwrite == True:
                print('Overwriting existing files.')
                self.clearDir(directory)
                self.directories.append(directory)
                return
        raise RuntimeError('Could not create directory: %s'%directory)

    def makeDirectories(self):
        '''
        Makes temporary directories for all processes.
        '''
        for procID in range(1, self.nproc + 1):
            directory = self.getProcDir(procID)
            self.makeDir(directory) 

    def makeInvIterDir(self):
        '''
        Makes directories for each iteration step for the output.
        '''
        invIterDir = self.cwd + '/it_%s'%(self.citer)
        err = os.system('mkdir %s'%invIterDir)
        if err is 256:
            if self.overwrite:
                self.clearDir(invIterDir)
        elif err is not 0:
            raise RuntimeError('Could not create directory: %s'%invIterDir)
        self.cInvIterDir = invIterDir        

    def clearDir(self, directory):
        '''
        Wipes a certain directory.
        '''
        #print('Wiping directory %s...'%directory)
        for filename in os.listdir(directory):
            filenp = os.path.join(directory, filename)
            os.remove(filenp)
        
    def clearDirectories(self):
        '''
        Wipes all generated temporary directories.
        '''
        for directory in self.directories:
            self.clearDir(directory) 

    def rmDir(self, directory):
        #print('Removing directory %s...'%directory)
        return os.rmdir(directory)

    def removeDirectories(self):
        '''
        Removes all generated temporary directories.
        '''
        for directory in self.directories:
            self.rmDir(directory) 
        self.directories = []
        
    def getProcDir(self, procID):
        '''
        Returns the temporary directory for a certain process
        with procID = process number.
        '''
        return os.path.join(self.cwd, self.folder) + str(procID)

    def getTraceIDs4Sources(self, sourceIDs):
        '''
        Returns corresponding trace IDs for a set of given source IDs.
        '''
        traceIDs = []
        for traceID in self.traces.keys():
            if self.traces[traceID]['source'] in sourceIDs:
                traceIDs.append(traceID)
        return traceIDs

    def getTraceIDs4Source(self, sourceID):
        '''
        Returns corresponding trace IDs for a source ID.
        '''
        traceIDs = []
        for traceID in self.traces.keys():
            if self.traces[traceID]['source'] == sourceID:
                traceIDs.append(traceID)
        return traceIDs

    def copyArrivals(self, target = None):
        '''
        Copies the FMM output file (self.ttim) to a specific target file.
        Default target is self.mtrav (model travel times).
        '''
        if target == None:
            target = os.path.join(self.cwd, self.mtrav)
        os.system('cp %s %s'%(os.path.join(
                    self.cInvIterDir, self.ttim), target))

    def saveVgrid(self):
        '''
        Saves the current velocity grid for the current iteration step.
        '''
        vgpath = os.path.join(self.cwd, self.cvg)
        os.system('cp %s %s'%(vgpath, self.cInvIterDir))

    def calcSrcPerKernel(self):
        '''
        (Equally) distributes all sources depending on the number of processes (kernels).
        Returns two integer values.
        First: minimum number of sources for each process
        Second: remaining sources (always less than number of processes)
        '''
        nsrc = self.readNsrc()
        if self.nproc > nsrc:
            print('Warning: Number of spawned processes higher than number of sources')
        return nsrc/self.nproc, nsrc%self.nproc

    def srcIDs4Kernel(self, procID):
        '''
        Calculates and returns all source IDs for a given process ID.
        '''
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
    
    def readNsrc(self):
        srcfile = open(self.csl, 'r')
        nsrc = int(srcfile.readline())
        srcfile.close()
        return nsrc

    def readNtraces(self):
        '''
        Reads the total number of traces from self.rec header.
        '''
        recfile = open(self.rec, 'r')
        nrec = int(recfile.readline())
        recfile.close()
        return nrec

    def readSrcFile(self):
        '''
        Reads the whole sourcefile and returns structured information in a dictionary.
        '''
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
        '''
        Reads the receiver input file and returns the information
        in a structured dictionary.
        '''
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

    def readArrivals(self, procID):
        '''
        Reads the arrival times from a temporary process directory,
        changes local to global sourceIDs and traceIDs and returns
        a list of arrival times.
        '''
        directory = self.getProcDir(procID)
        arrfile = open(os.path.join(directory, self.ttim), 'r')
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

    def readRays(self, procID):
        '''
        Reads rays output from a temporary process directory and returns
        the information in a structured dictionary.
        '''
        directory = self.getProcDir(procID)
        raysfile = open(directory + '/rays.dat', 'r')
        sourceIDs = self.srcIDs4Kernel(procID)
        
        rays = {}
        for sourceID in sourceIDs:
            traceIDs = self.getTraceIDs4Source(sourceID)
            for traceID in traceIDs:
                line1 = raysfile.readline().split()
                if line1 != []:
                    # recID and srcID for the individual processor will not be needed
                    recID_proc, srcID_proc, ray, normal, nsec = line1
                    raysecs = {}

                    for sec in range(int(nsec)):
                        line2 = raysfile.readline().split()
                        npoints, region, diff, head = line2
                        raypoints = []
                    
                        for j in range(int(npoints)):
                            raypoints.append(raysfile.readline())

                        raysecs[sec] = {'npoints': npoints,
                                        'region': region,
                                        'diff': diff,
                                        'head': head,
                                        'raypoints': raypoints
                                        }


                    rays[traceID] = {'sourceID': sourceID,
                                     'raypath': ray,
                                     'normal': normal,
                                     'nsec': nsec,
                                     'raysections': raysecs
                                     }
        return rays

    def writeSrcFile(self, procID):
        '''
        Writes a source input file for a process with ID = procID.
        '''
        directory = self.getProcDir(procID)
        srcfile = open(os.path.join(directory, self.csl), 'w')
        sourceIDs = self.srcIDs4Kernel(procID)

        srcfile.write('%s\n'%len(sourceIDs))
        for sourceID in sourceIDs:
            source = self.sources[sourceID]
            coords = source['coords']
            interactions = source['interactions']
            srcfile.write('%s\n'%source['teleflag'])
            srcfile.write('%s %s %s\n'%(float(coords[0]), float(coords[1]), float(coords[2])))
            srcfile.write('%s\n'%source['numpaths'])
            srcfile.write('%s\n'%source['steps'])
            srcfile.write('%s %s\n'%(int(interactions[0]), int(interactions[1])))
            srcfile.write('%s\n'%source['veltype'])

    def writeTracesFile(self, procID):
        '''
        Writes a receiver input file for a process with ID = procID.
        '''
        directory = self.getProcDir(procID)
        recfile = open('%s/receivers.in'%directory, 'w')
        sourceIDs = self.srcIDs4Kernel(procID)
        traceIDs = self.getTraceIDs4Sources(sourceIDs)

        recfile.write('%s\n'%len(traceIDs))
        for traceID in traceIDs:
            trace = self.traces[traceID]
            coords = trace['coords']
            source = int(trace['source']) - sourceIDs[0] + 1
            recfile.write('%s %s %s\n'%(float(coords[0]), float(coords[1]), float(coords[2])))
            recfile.write('%s\n'%trace['paths'])
            recfile.write('%s\n'%source)
            recfile.write('%s\n'%trace['path'])

    def mergeArrivals(self, directory):
        '''
        Merges the arrival times for all processes to self.cInvIterDir.
        '''
        arrfn = os.path.join(directory, self.ttim)
        arrivalsOut = open(arrfn, 'w')
        print('Merging %s...'%self.ttim)
        for procID in range(1, self.nproc + 1):
            arrivals = self.readArrivals(procID)
            for line in arrivals:
                arrivalsOut.write('%6s %6s %6s %6s %15s %5s %5s\n'%tuple(line))

        os.system('ln -fs %s %s'%(arrfn, os.path.join(self.cwd, self.ttim)))

    def mergeRays(self, directory):
        '''
        Merges the ray paths for all processes to self.cInvIterDir.
        '''
        print('Merging rays.dat...')
        with open(directory + '/rays.dat', 'w') as outfile:
            for procID in range(1, self.nproc + 1):
                rays = self.readRays(procID)
                for traceID in rays:
                    ray = rays[traceID]
                    outfile.write('%6s %6s %6s %6s %6s\n'%(traceID,
                                                           ray['sourceID'],
                                                           ray['raypath'],
                                                           ray['normal'],
                                                           ray['nsec']))
                    for sec in range(int(ray['nsec'])):
                        raysec = ray['raysections'][sec]
                        outfile.write('%6s %6s %6s %6s\n'%(raysec['npoints'],
                                                           raysec['region'],
                                                           raysec['diff'],
                                                           raysec['head']))
                        outfile.writelines(raysec['raypoints'])
                                                  
    def mergeFrechet(self, directory):
        '''
        Merges the frechet derivatives for all processes to self.cInvIterDir.
        '''
        frechfnout = os.path.join(directory, self.frechout)
        print('Merging %s...'%self.frechout)
        with open(frechfnout, 'w') as outfile:
            for procID in range(1, self.nproc + 1):
                filename = os.path.join(self.getProcDir(procID), self.frechout)
                with open(filename) as infile:
                    for sourceID in self.srcIDs4Kernel(procID):
                        for traceID in self.getTraceIDs4Source(sourceID):
                            recID_proc, srcID_proc, ray, normal, NPDEV = infile.readline().split()
                            outfile.write('%6s %6s %6s %6s %6s\n'%(traceID, sourceID, ray, normal, NPDEV))
                            for index in range(int(NPDEV)):
                                outfile.write(infile.readline())

        os.system('ln -fs %s %s'%(frechfnout, os.path.join(self.cwd, self.frechout)))

    def mergeOutput(self, directory):
        '''
        Calls self.mergeArrivals, self.mergeFrechet and self.mergeRays.
        '''
        self.mergeArrivals(directory)
        self.mergeFrechet(directory)
        self.mergeRays(directory)

    def unlink(self, filepath):
        return os.system('unlink %s'%filepath)
        
    def _printLine(self):
        print('----------------------------------------')

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
    outfile.write('# vtk DataFile Version 3.1\n')
    outfile.write('Velocity on FMTOMO vgrids.in points\n')
    outfile.write('ASCII\n')
    outfile.write('DATASET STRUCTURED_POINTS\n')

    outfile.write('DIMENSIONS %d %d %d\n' % (nX, nY, nZ))
    outfile.write('ORIGIN %f %f %f\n' % (sX, sY, sZ))
    outfile.write('SPACING %f %f %f\n' % (dX, dY, dZ))

    outfile.write('POINT_DATA %15d\n' % (nPoints))
    if absOrRel == 'abs':
        outfile.write('SCALARS velocity float %d\n' %(1))
    if absOrRel == 'relDepth':
        outfile.write('SCALARS velocity2depthMean float %d\n' %(1))
    elif absOrRel == 'rel':
        outfile.write('SCALARS velChangePercent float %d\n' % (1))
    outfile.write('LOOKUP_TABLE default\n')

    pointsPerR = nTheta * nPhi

    # write velocity
    if absOrRel == 'abs':
        print("Writing velocity values to VTK file...")
        for velocity in vel:
            outfile.write('%10f\n' %velocity)
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
                    outfile.write('%10f\n' %(vel - velmean))
                veldepth = []
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
            outfile.write('%10f\n' % velocity)
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
        if firstline == '': 
            break  # break at EOF
        fl_list = firstline.split()
        raynumber = int(fl_list[0])
        shotnumber = int(fl_list[1])
        rayValid = int(fl_list[4])  # is zero if the ray is invalid
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
        outfile.write('# vtk DataFile Version 3.1\n')
        outfile.write('FMTOMO rays\n')
        outfile.write('ASCII\n')
        outfile.write('DATASET POLYDATA\n')
        outfile.write('POINTS %15d float\n' % (nPoints))

        # write coordinates
        # print("Writing coordinates to VTK file...")
        for raynumber in rays[shotnumber].keys():
            for raypoint in rays[shotnumber][raynumber]:
                outfile.write('%10f %10f %10f \n' % (raypoint[0], raypoint[1], raypoint[2]))

        outfile.write('LINES %15d %15d\n' % (len(rays[shotnumber]), len(rays[shotnumber]) + nPoints))

        # write indices
        # print("Writing indices to VTK file...")
        count = 0
        for raynumber in rays[shotnumber].keys():
            outfile.write('%d ' % (len(rays[shotnumber][raynumber])))
            for index in range(len(rays[shotnumber][raynumber])):
                outfile.write('%d ' % (count))
                count += 1
            outfile.write('\n')


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
    outfile.write('%10s %10s \n' % (1, 1))
    outfile.write('%10s %10s %10s\n' % (nR, nTheta, nPhi))
    outfile.write('%10s %10s %10s\n' % (dR, np.deg2rad(dTheta), np.deg2rad(dPhi)))
    outfile.write('%10s %10s %10s\n' % (sR, np.deg2rad(sTheta), np.deg2rad(sPhi)))

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

                outfile.write('%10s %10s\n' % (velocity, decm))
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
    outfile.write('%10s %10s \n' % (1, 1))
    outfile.write('%10s %10s %10s\n' % (nR, nTheta, nPhi))
    outfile.write('%10s %10s %10s\n' % (dR, np.deg2rad(dTheta), np.deg2rad(dPhi)))
    outfile.write('%10s %10s %10s\n' % (sR, np.deg2rad(sTheta), np.deg2rad(sPhi)))

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

                outfile.write('%10s %10s\n' % (velocity, decm))
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
