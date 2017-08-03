#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import numpy
from pylot.core.pick.utils import getnoisewin

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--t', type=numpy.array, help='numpy array of time stamps')
    parser.add_argument('--t1', type=float, help='time from which relativ to it noise window is extracted')
    parser.add_argument('--tnoise', type=float, help='length of time window [s] for noise part extraction')
    parser.add_argument('--tgap', type=float, help='safety gap between signal (t1=onset) and noise')
    args = parser.parse_args()
    getnoisewin(args.t, args.t1, args.tnoise, args.tgap)
