#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obspy.core import read
from obspy.signal.trigger import coincidenceTrigger


class CoincidenceTimes(object):
    def __init__(self, st, comp='Z', coinum=4, sta=1., lta=10., on=5., off=1.):
        _type = 'recstalta'
        self.coinclist = self.createCoincTriggerlist(data=st, trigcomp=comp,
                                                     coinum=coinum, sta=sta,
                                                     lta=lta, trigon=on,
                                                     trigoff=off, type=_type)

    def __str__(self):
        n = 1
        out = ''
        for time in self.getCoincTimes():
            out += 'event no. {0}: starttime is {1}\n'.format(n, time)
            n += 1
        return out

    def getCoincTimes(self):
        timelist = []
        for info in self.getCoincList():
            timelist.append(info['time'])

        return timelist

    def getCoincList(self):
        return self.coinclist

    def createCoincTriggerlist(self, data, trigcomp, coinum, sta, lta,
                               trigon, trigoff, type):
        '''
        uses a coincidence trigger to detect all events in the given
        dataset
        '''

        triggerlist = coincidenceTrigger(type, trigon, trigoff,
                                         data.select(component=trigcomp),
                                         coinum, sta=sta, lta=lta)
        return triggerlist


def main():
    data = read('/data/SDS/2014/1A/ZV??/?H?.D/*.365')
    data.filter(type='bandpass', freqmin=5., freqmax=30.)
    coincs = CoincidenceTimes(data)
    print(coincs)


if __name__ == '__main__':
    main()
