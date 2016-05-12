#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
from obspy import UTCDateTime
import sys

def time_from_header(header):
    '''
    Function takes in the second line from a .gse file and takes out the date and time from that line.
    :param header: second line from .gse file
    :type header: string
    :return: a list of integers of form [year, month, day, hour, minute, second, microsecond]
    '''
    timeline = header.split(' ')
    time = timeline[1].split('/') + timeline[2].split(':')
    time = time[:-1] + time[-1].split('.')
    return [int(t) for t in time]

def check_time(datetime):
    '''
    Function takes in date and time as list and validates it's values by trying to make an UTCDateTime object from it
    :param datetime: list of integers [year, month, day, hour, minute, second, microsecond]
    :type datetime: list
    :return: returns True if Values are in supposed range, returns False otherwise
    '''
    try:
        UTCDateTime(*datetime)
        return True
    except ValueError:
        return False

def get_file_list(root_dir):
    '''
    Function uses a directorie to get all the *.gse files from it.
    :param root_dir: a directorie leading to the .gse files
    :type root_dir: string
    :return: returns a list of filenames (without path to them)
    '''
    file_list = glob.glob1(root_dir, '*.gse')
    return file_list

def checks_station_second(datetime):
    '''
    Function uses the given list to check if the parameter 'second' is set to 60 by mistake
     and sets the time correctly if so. Can only correct time if no date change would be necessary.
    :param datetime: [year, month, day, hour, minute, second, microsecond]
    :return: returns the input with the correct value for second
    '''
    if datetime[5] == 60:
        if datetime[4] == 59:
            if datetime[3] == 23:
                print 'Date should be next day.'+file
                raise ValueError
            else:
                datetime[3] += 1
                datetime[4] = 0
                datetime[5] = 0
        else:
            datetime[4] += 1
            datetime[5] = 0
    return datetime

def make_time_line(line, datetime):
    '''
    Function takes in the original line from a .gse file and a list of date and time values to make a new line with
    corrected date and time.
    :param line: second line from .gse file.
    :type line: string
    :param datetime: list of integers [year, month, day, hour, minute, second, microsecond]
    :type datetime: list
    :return: returns a string to write it into a file.
    '''
    insertion = '{:02d}'.format(int(datetime[3])) + ':' + '{:02d}'.format(int(datetime[4])) + ':' + '{:02d}'.format(int(datetime[5])) + '.000'
    newline = line[:16]+insertion+line[28:]
    return newline

def evt_head_check(root_dir,out_dir = None):
    '''
    A function to make sure that an arbitrary number of .gse files have correct values in their header.
    :param root_dir: a directory leading to the .gse files.
    :type root_dir: string
    :param out_dir: a directory to store the new files somwhere els.
    :return: returns nothing
    '''
    if not out_dir:
        print 'WARNING files are going to be overwritten!'
        inp = str(raw_input('Continue? [y/n]'))
        if inp == 'y':
            pass
        else:
            sys.exit()
    Filelist = get_file_list(root_dir)
    debugcounter = 0
    for i in range(len(Filelist)):
        inFile = open(root_dir+'/'+Filelist[i], 'r')
        lines = inFile.readlines()
        datetime = time_from_header(lines[1])
        if check_time(datetime):
            continue
        else:
            debugcounter += 1
            datetime = checks_station_second(datetime)
            print 'writing ' + Filelist[i]
            timeline = make_time_line(lines[1],datetime)
            # write File
            lines[1] = timeline
            if not out_dir:
                out = open(root_dir+Filelist[i], 'w')
                out.writelines(lines)
                out.close()
            else:
                out = open(out_dir+Filelist[i], 'w')
                out.writelines(lines)
                out.close()
        inFile.close()
    print debugcounter