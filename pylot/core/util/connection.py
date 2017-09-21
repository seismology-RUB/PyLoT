#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from urllib2 import urlopen
except:
    from urllib.request import urlopen


def checkurl(url='https://ariadne.geophysik.ruhr-uni-bochum.de/trac/PyLoT/'):
    try:
        urlopen(url, timeout=1)
        return True
    except:
        pass
    return False
