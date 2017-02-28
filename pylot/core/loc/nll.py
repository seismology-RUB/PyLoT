#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import os
import glob
from obspy import read_events
from pylot.core.io.phases import writephases
from pylot.core.util.utils import getPatternLine, runProgram, which
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()

class NLLocError(EnvironmentError):
    pass

def export(picks, fnout, parameter):
    '''
    Take <picks> dictionary and exports picking data to a NLLOC-obs
    <phasefile> without creating an ObsPy event object.

    :param picks: picking data dictionary
    :type picks: dict

    :param fnout: complete path to the exporting obs file
    :type fnout: str
 
    :param: parameter, all input information
    :type:  object
    '''
    # write phases to NLLoc-phase file
    writephases(picks, 'NLLoc', fnout, parameter)


def modify_inputs(ctrfn, root, nllocoutn, phasefn, tttn):
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


def locate(fnin):
    """
    takes an external program name
    :param fnin:
    :return:
    """

    exe_path = which('NLLoc')
    if exe_path is None:
        raise NLLocError('NonLinLoc executable not found; check your '
                         'environment variables')

    # locate the event utilizing external NonLinLoc installation
    try:
        runProgram(exe_path, fnin)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(e.output)


def read_location(fn):
    path, file = os.path.split(fn)
    file = glob.glob1(path, file +  '.[0-9]*.grid0.loc.hyp')
    if len(file) > 1:
        raise IOError('ambiguous location name {0}'.format(file))
    fn = os.path.join(path, file[0])
    return read_events(fn)[0]


if __name__ == '__main__':
    pass
