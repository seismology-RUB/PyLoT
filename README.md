# PyLoT

version: 0.3

The Python picking and Localisation Tool

This python library contains a graphical user interfaces for picking seismic phases. This software needs [ObsPy][ObsPy]
and the PySide2 Qt5 bindings for python to be installed first.

PILOT has originally been developed in Mathworks' MatLab. In order to distribute PILOT without facing portability
problems, it has been decided to redevelop the software package in Python. The great work of the ObsPy group allows easy
handling of a bunch of seismic data and PyLoT will benefit a lot compared to the former MatLab version.

The development of PyLoT is part of the joint research project MAGS2 and AlpArray.

## Installation

At the moment there is no automatic installation procedure available for PyLoT. Best way to install is to clone the
repository and add the path to your Python path.

It is highly recommended to use Anaconda for a simple creation of a Python installation using either the *pylot.yml* or the *requirements.txt* file found in the PyLoT root directory. First make sure that the *conda-forge* channel is available in your Anaconda installation:

    conda config --add channels conda-forge

Afterwards run (from the PyLoT main directory)

    conda create --name pylot_38 --file requirements.txt
or
    
    conda env create -f pylot.yml

#### Prerequisites:

In order to run PyLoT you need to install:

- Python 3
- obspy
- pyside2
- pyqtgraph
- cartopy

(the following are already dependencies of the above packages):
- scipy
- numpy
- matplotlib

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

and some extra information on error estimates (just needed for reading old PILOT data) and the Richter magnitude scaling
relation

    cp path-to-pylot/inputs/PILOT_TimeErrors.in path-to-pylot/inputs/richter_scaling.data ~/.pylot/

You may need to do some modifications to these files. Especially folder names should be reviewed.

PyLoT has been tested on Mac OSX (10.11), Debian Linux 8 and on Windows 10.

## Release notes

#### Features:

- event organisation in project files and waveform visualisation
- consistent manual phase picking through predefined SNR dependant zoom level
- consistent automatic phase picking routines using Higher Order Statistics, AIC and Autoregression
- interactive tuning of auto-pick parameters
- uniform uncertainty estimation from waveform's properties for automatic and manual picks
- pdf representation and comparison of picks taking the uncertainty intrinsically into account
- Richter and moment magnitude estimation
- location determination with external installation of [NonLinLoc](http://alomax.free.fr/nlloc/index.html)

#### Known issues:

We hope to solve these with the next release.

## Staff

Original author(s): L. Kueperkoch, S. Wehling-Benatelli, M. Bischoff (PILOT)

Developer(s): S. Wehling-Benatelli, L. Kueperkoch, K. Olbert, M. Bischoff, C. Wollin, M. Rische, M. Paffrath

Others: A. Bruestle, T. Meier, W. Friederich


[ObsPy]: http://github.com/obspy/obspy/wiki

September 2017
