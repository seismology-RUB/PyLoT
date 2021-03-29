#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

try:
    import pyqtgraph as pg
except Exception as e:
    print('Warning: Could not import module pyqtgraph.')
try:
    from PySide import QtCore
except Exception as e:
    print('Warning: Could not import module QtCore.')


from pylot.core.util.utils import pick_color


def pick_linestyle_pg(picktype, key):
    """
    Get Qt line style by  picktype and pick parameter (earliest/latest possible pick, symmetric picking error or
    most probable pick)
    :param picktype: 'manual' or 'automatic'
    :type picktype: str
    :param key: which pick parameter should be plotted, 'mpp', 'epp', 'lpp' or 'spe'
    :type key: str
    :return: Qt line style parameters
    :rtype:
    """
    linestyles_manu = {'mpp': (QtCore.Qt.SolidLine, 2.),
                       'epp': (QtCore.Qt.DashLine, 1.),
                       'lpp': (QtCore.Qt.DashLine, 1.),
                       'spe': (QtCore.Qt.DashLine, 1.)}
    linestyles_auto = {'mpp': (QtCore.Qt.DotLine, 2.),
                       'epp': (QtCore.Qt.DashDotLine, 1.),
                       'lpp': (QtCore.Qt.DashDotLine, 1.),
                       'spe': (QtCore.Qt.DashDotLine, 1.)}
    linestyles = {'manual': linestyles_manu,
                  'auto': linestyles_auto}
    return linestyles[picktype][key]


def which(program, parameter):
    """
    takes a program name and returns the full path to the executable or None
    modified after: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    :param program: name of the desired external program
    :type program: str
    :return: full path of the executable file
    :rtype: str
    """
    try:
        from PySide.QtCore import QSettings
        settings = QSettings()
        for key in settings.allKeys():
            if 'binPath' in key:
                os.environ['PATH'] += ':{0}'.format(settings.value(key))	
        nllocpath = ":" + parameter.get('nllocbin')
        os.environ['PATH'] += nllocpath
    except Exception as e:
        print(e.message)

    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program) 
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate

    return None


def make_pen(picktype, phase, key, quality):
    """
    Make  PyQtGraph.QPen
    :param picktype: 'manual' or 'automatic'
    :type picktype: str
    :param phase: 'P' or 'S'
    :type phase: str
    :param key: 'mpp', 'epp', 'lpp' or 'spe', (earliest/latest possible pick, symmetric picking error or
     most probable pick)
    :type key: str
    :param quality: quality class of pick, decides color modifier
    :type quality: int
    :return: PyQtGraph QPen
    :rtype: `~QPen`
    """
    if pg:
        rgba = pick_color(picktype, phase, quality)
        linestyle, width = pick_linestyle_pg(picktype, key)
        pen = pg.mkPen(rgba, width=width, style=linestyle)
        return pen

