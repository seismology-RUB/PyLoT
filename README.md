# PyLoT

version: 0.3

The Python picking and Localisation Tool

This python library contains a graphical user interfaces for picking seismic phases. This software needs [ObsPy][ObsPy]
and the PySide2 Qt5 bindings for python to be installed first.

PILOT has originally been developed in Mathworks' MatLab. In order to distribute PILOT without facing portability
problems, it has been decided to redevelop the software package in Python. The great work of the ObsPy group allows easy
handling of a bunch of seismic data and PyLoT will benefit a lot compared to the former MatLab version.

The development of PyLoT is part of the joint research project MAGS2, AlpArray and AdriaArray.

## Installation

At the moment there is no automatic installation procedure available for PyLoT. Best way to install is to clone the
repository and add the path to your Python path.

It is highly recommended to use Anaconda for a simple creation of a Python installation using either the *pylot.yml* or the *requirements.txt* file found in the PyLoT root directory. First make sure that the *conda-forge* channel is available in your Anaconda installation:

    conda config --add channels conda-forge

Afterwards run (from the PyLoT main directory where the files *requirements.txt* and *pylot.yml* are located)

    conda env create -f pylot.yml
or
    
    conda create -c conda-forge --name pylot_311 python=3.11 --file requirements.txt

to create a new Anaconda environment called *pylot_311*.

Afterwards activate the environment by typing

    conda activate pylot_311

#### Prerequisites:

In order to run PyLoT you need to install:

- Python 3
- cartopy
- joblib
- obspy
- pyaml
- pyqtgraph
- pyside2

(the following are already dependencies of the above packages):
- scipy
- numpy
- matplotlib

#### Some handwork:

Some extra information on error estimates (just needed for reading old PILOT data) and the Richter magnitude scaling
relation

    cp path-to-pylot/inputs/PILOT_TimeErrors.in path-to-pylot/inputs/richter_scaling.data ~/.pylot/

You may need to do some modifications to these files. Especially folder names should be reviewed.

PyLoT has been tested on Mac OSX (10.11), Debian Linux 8 and on Windows 10/11.

## Example Dataset
An example dataset with waveform data, metadata and automatic picks in the obspy-dmt dataset format for testing the teleseismic picking can be found at https://zenodo.org/doi/10.5281/zenodo.13759803

## Release notes

#### Features:

- event organisation in project files and waveform visualisation
- consistent manual phase picking through predefined SNR dependant zoom level
- consistent automatic phase picking routines using Higher Order Statistics, AIC and Autoregression
- pick correlation correction for teleseismic waveforms
- interactive tuning of auto-pick parameters
- uniform uncertainty estimation from waveform's properties for automatic and manual picks
- pdf representation and comparison of picks taking the uncertainty intrinsically into account
- Richter and moment magnitude estimation
- location determination with external installation of [NonLinLoc](http://alomax.free.fr/nlloc/index.html)

#### Known issues:

Current release is still in development progress and has several issues. We are currently lacking manpower, but hope to assess many of the issues in the near future.

## Staff

Original author(s): M. Rische, S. Wehling-Benatelli, L. Kueperkoch, M. Bischoff (PILOT)

Developer(s): S. Wehling-Benatelli, M. Paffrath, L. Kueperkoch, K. Olbert, M. Bischoff, C. Wollin, M. Rische, D. Arnold, K. Cökerim, S. Zimmermann

Others: A. Bruestle, T. Meier, W. Friederich


[ObsPy]: http://github.com/obspy/obspy/wiki

September 2024