from distutils.core import setup

setup(
    name='PyLoT',
    version='0.1a1',
    packages=['pylot', 'pylot.core', 'pylot.core.loc', 'pylot.core.pick',
              'pylot.core.read', 'pylot.core.util', 'pylot.core.active',
              'pylot.core.analysis', 'pylot.testing'],
    url='dummy',
    license='LGPLv3',
    author='Sebastian Wehling-Benatelli',
    author_email='sebastian.wehling@rub.de',
    description='Comprehensive Python picking and Location Toolbox for seismological data.'
)
