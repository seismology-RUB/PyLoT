#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pylot.core.io.phases import reassess_pilot_db
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'S. Wehling-Benatelli'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='reassess old PILOT event data base in terms of consistent '
                    'automatic uncertainty estimation',
        epilog='Script written by {author} belonging to PyLoT version'
               ' {version}\n'.format(author=__author__,
                                     version=__version__)
    )

    parser.add_argument(
        'root', type=str, help='specifies the root directory'
    )
    parser.add_argument(
        'db', type=str, help='specifies the database name'
    )
    parser.add_argument(
        '--output', '-o', type=str, help='path to the output directory',
        dest='output'
    )
    parser.add_argument(
        '--parameterfile', '-p', type=str,
        help='full path to the parameterfile', dest='parfile'
    )
    parser.add_argument(
        '--verbosity', '-v', action='count', help='increase output verbosity',
        default=0, dest='verbosity'
    )

    args = parser.parse_args()
    reassess_pilot_db(args.root, args.db, args.output, args.parfile, args.verbosity)
