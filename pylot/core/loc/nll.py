#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import os
from obspy.core.event import readEvents
from pylot.core.pick.utils import writephases
from pylot.core.util.utils import getPatternLine
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()

def picksExport(picks, phasefile):
    '''
    Take <picks> dictionary and exports picking data to a NLLOC-obs <phasefile> without creating an ObsPy event object.
    :param picks: picking data dictionary
    :type picks: dict
    :param phasefile: complete path to the exporting obs file
    :type phasefile: str
    '''
    # write phases to NLLoc-phase file
    writephases(picks, 'NLLoc', phasefile)

def modfiyInputFile(fn, root, outpath, phasefn, tttn):
    '''

    :param fn:
    :param root:
    :param outpath:
    :param phasefile:
    :param tttable:
    :return:
    '''
    # For locating the event we have to modify the NLLoc-control file!
    # create comment line for NLLoc-control file NLLoc-output file
    print ("Modifying  NLLoc-control file %s ..." % fn)
    nllocout = os.path.join(root,'loc', outpath)
    phasefile = os.path.join(root, 'obs', phasefn)
    tttable = os.path.join(root, 'time', tttn)
    locfiles = 'LOCFILES %s NLLOC_OBS %s %s 0\n' % (phasefile, tttable, nllocout)

    # modification of NLLoc-control file
    curlocfiles = getPatternLine(fn, 'LOCFILES')
    nllfile = open(fn, 'r')
    filedata = nllfile.read()
    if filedata.find(locfiles) < 0:
        # replace old command
        filedata = filedata.replace(curlocfiles, locfiles)
        nllfile = open(fn, 'w')
        nllfile.write(filedata)
    nllfile.close()

def locate(call, fnin):
    '''
    Takes paths to NLLoc executable <call> and input parameter file <fnin> and starts the location calculation.
    :param call: full path to NLLoc executable
    :type call: str
    :param fnin: full path to input parameter file
    :type fnin: str
    '''
    # locate the event
    subprocess.call([call, fnin])

def readLocation(fn):
    pass

if __name__=='__main__':
    pass
