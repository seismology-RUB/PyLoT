#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from pylot.core.util.version import get_git_version as _getVersionString
from pylot.core.io.phases import reasses_pilot_event

__version__ = _getVersionString()
__author__ = 'sebastianw'

def reassess_pilot_event():
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    args = parser.parse_args()
    reasses_pilot_event(args.id)
