#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import os
from pylot.core.io.phases import writephases
from pylot.core.util.utils import getPatternLine
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()


def picksExport(picks, locrt, phasefile):
    '''
    Take <picks> dictionary and exports picking data to a NLLOC-obs
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param locrt: choose location routine
    :type locrt: str

    :param phasefile: complete path to the exporting obs file
    :type phasefile: str
    '''
    # write phases to NLLoc-phase file
    writephases(picks, locrt, phasefile)


def modifyInputFile(ctrfn, root, nllocoutn, phasefn, tttn):
    '''
    :param ctrfn: name of NLLoc-control file
    :type: str

    :param root: root path to NLLoc working directory
    :type: str

    :param nllocoutn: name of NLLoc-location output file
    :type: str

    :param phasefn: name of NLLoc-input phase file
    :type: str

    :param tttn: pattern of precalculated NLLoc traveltime tables
    :type: str
    '''
    # For locating the event the NLLoc-control file has to be modified!
    # create comment line for NLLoc-control file NLLoc-output file
    ctrfile = os.path.join(root, 'run', ctrfn)
    nllocout = os.path.join(root, 'loc', nllocoutn)
    phasefile = os.path.join(root, 'obs', phasefn)
    tttable = os.path.join(root, 'time', tttn)
    locfiles = 'LOCFILES %s NLLOC_OBS %s %s 0\n' % (phasefile, tttable, nllocout)

    # modification of NLLoc-control file
    print ("Modifying  NLLoc-control file %s ..." % ctrfile)
    curlocfiles = getPatternLine(ctrfile, 'LOCFILES')
    nllfile = open(ctrfile, 'r')
    filedata = nllfile.read()
    if filedata.find(locfiles) < 0:
        # replace old command
        filedata = filedata.replace(curlocfiles, locfiles)
        nllfile = open(ctrfile, 'w')
        nllfile.write(filedata)
    nllfile.close()


def locate(call, fnin):
    '''
    Takes paths to NLLoc executable <call> and input parameter file <fnin>
    and starts the location calculation.

    :param call: full path to NLLoc executable
    :type call: str

    :param fnin: full path to input parameter file
    :type fnin: str
    '''

    # locate the event
    subprocess.call([call, fnin])


def readLocation(fn):
    pass


if __name__ == '__main__':
    pass
