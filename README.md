# PyLoT

version: 0.2

The Python picking and Localisation Tool

This python library contains a graphical user interfaces for picking
seismic phases. This software needs [ObsPy][ObsPy]
and the PySide Qt4 bindings for python to be installed first.

PILOT has originally been developed in Mathworks' MatLab. In order to
distribute PILOT without facing portability problems, it has been decided
to redevelop the software package in Python. The great work of the ObsPy
group allows easy handling of a bunch of seismic data and PyLoT will
benefit a lot compared to the former MatLab version.

The development of PyLoT is part of the joint research project MAGS2 and AlpArray.

## Installation

At the moment there is no automatic installation procedure available for PyLoT.
Best way to install is to clone the repository and add the path to your Python path.

#### Prerequisites:

In order to run PyLoT you need to install:

- python 2 or 3
- scipy
- numpy
- matplotlib
- obspy
- pyside

#### Some handwork:

PyLoT needs a properties folder on your system to work. It should be situated in your home directory 
(on Windows usually C:/Users/*username*):

    mkdir ~/.pylot

In the next step you have to copy some files to this directory:

*for local distance seismicity*

    cp path-to-pylot/inputs/pylot_local.in ~/.pylot/pylot.in

*for regional distance seismicity*

    cp path-to-pylot/inputs/pylot_regional.in ~/.pylot/pylot.in

*for global distance seismicity*

    cp path-to-pylot/inputs/pylot_global.in ~/.pylot/pylot.in

and some extra information on error estimates (just needed for reading old PILOT data) and the Richter magnitude scaling relation

    cp path-to-pylot/inputs/PILOT_TimeErrors.in path-to-pylot/inputs/richter_scaling.data ~/.pylot/

You may need to do some modifications to these files. Especially folder names should be reviewed.

PyLoT has been tested on Mac OSX (10.11), Debian Linux 8 and on Windows 10.


## Release notes

#### Features:

- centralize all functionalities of PyLoT and control them from within the main GUI
- handling multiple events inside GUI with project files (save and load work progress)
- GUI based adjustments of pick parameters and I/O
- interactive tuning of parameters from within the GUI
- call automatic picking algorithm from within the GUI
- comparison of automatic with manual picks for multiple events using clear differentiation of manual picks into 'tune' and 'test-set' (beta)
- manual picking of different (user defined) phase types
- phase onset estimation with ObsPy TauPy
- interactive zoom/scale functionalities in all plots (mousewheel, pan, pan-zoom)
- array map to visualize stations and control onsets (beta feature, switch to manual picks not implemented)

##### Platform support:
- Python 3 support
- Windows support

##### Performance:
- multiprocessing for automatic picking and restitution of multiple stations
- use pyqtgraph library for better performance on main waveform plot

##### Visualization:
- pick uncertainty (quality classes) visualization with gradients
- pick color unification for all plots
- new icons and stylesheets

#### Known Issues:
- some Qt related errors might occur at runtime
- filter toggle not working in pickDlg
- PyLoT data structure requires at least three parent directories for waveform data directory

## Staff

Original author(s): L. Kueperkoch, S. Wehling-Benatelli, M. Bischoff (PILOT)

Developer(s): S. Wehling-Benatelli, L. Kueperkoch, M. Paffrath, K. Olbert,
M. Bischoff, C. Wollin, M. Rische

Others: A. Bruestle, T. Meier, W. Friederich


[ObsPy]: http://github.com/obspy/obspy/wiki

September 2017
