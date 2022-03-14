#!/usr/bin/env python
# -*- coding: utf-8 -*-
from distutils.core import setup

setup(
    name='PyLoT',
    version='0.2',
    packages=['pylot', 'pylot.core', 'pylot.core.io', 'pylot.core.loc', 'pylot.core.pick', 'pylot.core.util',
              'pylot.core.analysis', 'pylot.styles'],
    requires=['obspy', 'PySide2', 'matplotlib', 'numpy', 'scipy', 'pyqtgraph', 'cartopy'],
    url='dummy',
    license='LGPLv3',
    author='Sebastian Wehling-Benatelli',
    author_email='sebastian.wehling@rub.de',
    description='Comprehensive Python picking and Location Toolbox for seismological data.'
)
