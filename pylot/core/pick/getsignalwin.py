#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy
from pylot.core.pick.utils import getsignalwin

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--t', type=~numpy.array, help='numpy array of time stamps')
    parser.add_argument('--t1', type=float, help='time from which relativ to it signal window is extracted')
    parser.add_argument('--tsignal', type=float, help='length of time window [s] for signal part extraction')
    args = parser.parse_args()
    getsignalwin(args.t, args.t1, args.tsignal)
