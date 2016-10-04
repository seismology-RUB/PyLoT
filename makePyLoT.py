#!/usr/bin/env python
# encoding: utf-8
from __future__ import print_function

"""
makePyLoT -- build and install PyLoT

makePyLoT is a python make file in order to establish the folder structure and
meet requisites

It defines
:class CLIError:
:method main:

:author:     Sebastian Wehling-Benatelli

:copyright:  2014 MAGS2 EP3 Working Group. All rights reserved.

:license:    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)

:contact:    sebastian.wehling@rub.de

updated: Updated
"""

import glob
import os
import sys
import shutil
import copy

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2014-11-26'
__updated__ = '2016-04-28'

DEBUG = 0
TESTRUN = 0
PROFILE = 0


class CLIError(Exception):
    """Generic exception to raise and log different fatal errors."""

    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg


def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = 'makePyLoT %s (%s)' % (
        program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''{0:s}

    Created by Sebastian Wehling-Benatelli on {1:s}.
    Copyright 2014 MAGS2 EP3 Working Group. All rights reserved.

    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)

    Distributed on an "AS IS" basis without warranties
    or conditions of any kind, either express or implied.

    USAGE
    '''.format(program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license,
                                formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-b", "--build", dest="build", action="store_true",
                            help="build PyLoT")
        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                            help="set verbosity level")
        parser.add_argument("-i", "--install", dest="install",
                            action="store_true",
                            help="install PyLoT on the system")
        parser.add_argument("-d", "--directory", dest="directory",
                            help="installation directory", metavar="RE")
        parser.add_argument('-V', '--version', action='version',
                            version=program_version_message)

        # Process arguments
        args = parser.parse_args()

        verbose = args.verbose
        build = args.build
        install = args.install
        directory = args.directory

        if verbose > 0:
            print("Verbose mode on")
        if install and not directory:
            raise CLIError("""Trying to install without appropriate
                           destination; please specify an installation
                           directory!""")
        if build and install:
            print("Building and installing PyLoT ...\n")
            buildPyLoT(verbose)
            installPyLoT(verbose)
        elif build and not install:
            print("Building PyLoT without installing! Please wait ...\n")
            buildPyLoT(verbose)
        cleanUp()
        return 0
    except KeyboardInterrupt:
        cleanUp(1)
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise e
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


def buildPyLoT(verbosity=None):
    system = sys.platform
    if verbosity > 1:
        msg = ("... on system: {0}\n"
               "\n"
               "              Current working directory: {1}\n"
               ).format(system, os.getcwd())
        print(msg)
    if system.startswith(('win', 'microsoft')):
        raise CLIError(
            "building on Windows system not tested yet; implementation pending")
    elif system == 'darwin':
        # create a symbolic link to the desired python interpreter in order to
        # display the right application name
        for path in os.getenv('PATH').split(':'):
            found = glob.glob(os.path.join(path, 'python'))
            if found:
                os.symlink(found, './PyLoT')
                break


def installPyLoT(verbosity=None):
    files_to_copy = {'autoPyLoT_local.in':['~', '.pylot'],
                     'autoPyLoT_regional.in':['~', '.pylot'],
                     'filter.in':['~', '.pylot']}
    if verbosity > 0:
        print ('starting installation of PyLoT ...')
    if verbosity > 1:
        print ('copying input files into destination folder ...')
    ans = input('please specify scope of interest '
                '([0]=local, 1=regional) :') or 0
    if not isinstance(ans, int):
        ans = int(ans)
    ans = 'local' if ans is 0 else 'regional'
    link_dest = []
    for file, destination in files_to_copy.items():
        link_file = ans in file
        if link_file:
            link_dest = copy.deepcopy(destination)
            link_dest.append('autoPyLoT.in')
            link_dest = os.path.join(*link_dest)
        destination.append(file)
        destination = os.path.join(*destination)
        srcfile = os.path.join('input', file)
        assert not os.path.isabs(srcfile), 'source files seem to be ' \
                                           'corrupted ...'
        if verbosity > 1:
            print ('copying file {file} to folder {dest}'.format(file=file, dest=destination))
        shutil.copyfile(srcfile, destination)
        if link_file:
            if verbosity:
                print('linking input file for autoPyLoT ...')
            os.symlink(destination, link_dest)




def cleanUp(verbosity=None):
    if verbosity >= 1:
        print('cleaning up build files...')
    if sys.platform == 'darwin':
        os.remove('./PyLoT')


if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
    if TESTRUN:
        import doctest

        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats

        profile_filename = 'makePyLoT_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
