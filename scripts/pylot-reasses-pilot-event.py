#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pylot.core.util.version import get_git_version as _getVersionString
from pylot.core.io.phases import reassess_pilot_event

__version__ = _getVersionString()
__author__ = 'sebastianw'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--directory', '-d', type=str, help='specifies the root directory (in '
                                            'most cases PILOT database folder)'
    )
    parser.add_argument(
        '--eventid', '-i', type=str, help='PILOT event identifier'
    )

    args = parser.parse_args()
    reassess_pilot_event(args.dir, args.id)
