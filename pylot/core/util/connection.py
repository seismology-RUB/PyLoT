#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    # noinspection PyUnresolvedReferences
    from urllib2 import urlopen
except:
    from urllib.request import urlopen


def checkurl(url='https://git.geophysik.ruhr-uni-bochum.de/marcel/pylot/'):
    """
    check if URL is available
    :param url: url
    :type url: str
    :return: available: True/False
    :rtype: bool
    """
    try:
        urlopen(url, timeout=1)
        return True
    except:
        pass
    return False
