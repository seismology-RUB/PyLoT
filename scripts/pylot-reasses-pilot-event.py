#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from pylot.core.pick.utils import reassess_pilot_event
from pylot.core.util.version import get_git_version as _getVersionString

__version__ = _getVersionString()
__author__ = 'sebastianw'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--TSNR', type=tuple, help='length of time windows around pick used to determine SNR \
        [s] (Tnoise, Tgap, Tsignal)', action=store_value
    )

    args = parser.parse_args()
    reassess_pilot_event(args.db, args.id, args.TSNR)
