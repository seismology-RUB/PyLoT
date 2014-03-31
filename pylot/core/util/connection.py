# -*- coding: utf-8 -*-

import urllib2


def checkurl(url='https://ariadne.geophysik.rub.de/trac/PyLoT'):
    try:
        urllib2.urlopen(url, timeout=1)
        return True
    except urllib2.URLError:
        pass
    return False
