#!/usr/bin/env python
# -*- coding: utf-8 -*-
from distutils.core import setup

setup(
    name='PyLoT',
    version='0.1a1',
    packages=['pylot', 'pylot.core', 'pylot.core.loc', 'pylot.core.pick',
              'pylot.core.io', 'pylot.core.util', 'pylot.core.active',
              'pylot.core.analysis', 'pylot.testing'],
    requires=['obspy', 'PySide', 'matplotlib', 'numpy'],
    url='dummy',
    license='LGPLv3',
    author='Sebastian Wehling-Benatelli',
    author_email='sebastian.wehling@rub.de',
    description='Comprehensive Python picking and Location Toolbox for seismological data.'
)
