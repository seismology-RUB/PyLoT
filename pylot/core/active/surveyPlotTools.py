# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
plt.interactive(True)

class regions(object):
    def __init__(self, ax, shot_dict):
        self.ax = ax
        self.shot_dict = shot_dict
        self._x0 = []
        self._y0 = []
        self._x1 = []
        self._y1 = []
        self.shots_found = {}
        self.shots_for_deletion = {}

    def _onselect(self, eclick, erelease):
        'eclick and erelease are matplotlib events at press and release'                                                                                            #print ' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata)
        #print ' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata)
        print 'region selected x0, y0 = (%3s, %3s), x1, y1 = (%3s, %3s)'%(eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata)
        x0 = min(eclick.xdata, erelease.xdata)
        x1 = max(eclick.xdata, erelease.xdata)
        y0 = min(eclick.ydata, erelease.ydata)
        y1 = max(eclick.ydata, erelease.ydata)
        self._x0.append(x0)
        self._x1.append(x1)
        self._y0.append(y0)
        self._y1.append(y1)
        self.markCurrentRegion(x0, x1, y0, y1)

    def chooseRectangles(self):
        from matplotlib.widgets import RectangleSelector

        print 'Select rectangle is active'
        return RectangleSelector(self.ax, self._onselect)

    def _getx0(self):
        return self._x0

    def _getx1(self):
        return self._x1

    def _gety0(self):
        return self._y0

    def _gety1(self):
        return self._y1

    def getShotDict(self):
        return self.shot_dict

    def getShotsForDeletion(self):
        return self.shots_for_deletion

    def findTracesInShotDict(self, picks = 'normal'):
        '''
        Returns traces corresponding to a certain area in a plot with all picks over the distances.
        '''
        print "findTracesInShotDict: Searching for marked traces in the shot dictionary... "

        for shot in self.shot_dict.values():
            whichpicks = {'normal': shot.getPick,
                          'includeCutOut': shot.getPick_backup}
            for index in range(len(self._getx1())):
                distancebin = (self._getx0()[index], self._getx1()[index])
                pickbin = (self._gety0()[index], self._gety1()[index])
                if shot.getTraceIDs4Dist(distancebin = distancebin) is not None:
                    for traceID in shot.getTraceIDs4Dist(distancebin = distancebin):
                        if pickbin[0] < whichpicks[picks](traceID) < pickbin[1]:
                            self.highlightPick(shot, traceID)
                            if shot.getShotnumber() not in self.shots_found.keys():
                                self.shots_found[shot.getShotnumber()] = []
                            if traceID not in self.shots_found[shot.getShotnumber()]:
                                self.shots_found[shot.getShotnumber()].append(traceID)
        self.refreshFigure()
        print self.shots_found

    def highlightPick(self, shot, traceID, annotations = True):
        self.ax.scatter(shot.getDistance(traceID), shot.getPick(traceID), s = 50, marker = 'o', facecolors = 'none', edgecolors = 'm', alpha = 1)
        if annotations == True:
            self.ax.annotate(s = 's%s|t%s'%(shot.getShotnumber(), traceID), xy = (shot.getDistance(traceID), shot.getPick(traceID)), fontsize = 'xx-small')
        self.ax.set_ylim(shot.getCut())

    def plotTracesInRegion(self):
        count = 0
        maxfigures = 20
        # if len(self.shots_found) == 0:
        self.findTracesInShotDict()

        if len(self.shots_found) > 0:
            for shot in self.shot_dict.values():
                for shotnumber in self.shots_found:
                    if shot.getShotnumber() == shotnumber:
                        for traceID in self.shots_found[shotnumber]:
                            count += 1
                            if count > maxfigures:
                                print 'Maximum number of figures (%s) reached. %sth figure was not opened.' %(maxfigures, count)
                                break
                            shot.plot_traces(traceID)
        else:
            print 'No picks yet defined in the regions x = (%s, %s), y = (%s, %s)' %(self._x0, self._x1, self._y0, self._y1)

    def plotTracesInRegion_withCutOutTraces(self):
        count = 0
        maxfigures = 20
        # if len(self.shots_found) == 0:
        self.findTracesInShotDict(picks = 'includeCutOut')

        if len(self.shots_found) > 0:
            for shot in self.shot_dict.values():
                for shotnumber in self.shots_found:
                    if shot.getShotnumber() == shotnumber:
                        for traceID in self.shots_found[shotnumber]:
                            count += 1
                            if count > maxfigures:
                                print 'Maximum number of figures (%s) reached. %sth figure was not opened.' %(maxfigures, count)
                                break
                            shot.plot_traces(traceID)
        else:
            print 'No picks yet defined in the regions x = (%s, %s), y = (%s, %s)' %(self._x0, self._x1, self._y0, self._y1)


    def setCurrentRegionsForDeletion(self):
        # if len(self.shots_found) == 0:
        self.findTracesInShotDict()

        for shotnumber in self.shots_found:
            if not shotnumber in self.shots_for_deletion:
                self.shots_for_deletion[shotnumber] = []
            for traceID in self.shots_found[shotnumber]:
                if not traceID in self.shots_for_deletion[shotnumber]:
                    self.shots_for_deletion[shotnumber].append(traceID)
        self.markAllRegions(color = 'red')
        print 'Marked regions for deletion'

    def markAllRegions(self, color = 'grey'):
        from matplotlib.patches import Rectangle

        for index in range(len(self._getx0())):
            x0 = self._getx0()[index]
            y0 = self._gety0()[index]
            x1 = self._getx1()[index]
            y1 = self._gety1()[index]

            self.ax.add_patch(Rectangle((x0, y0), (x1 - x0), (y1 - y0), alpha=0.5, facecolor = color))
            self.refreshFigure()

    def markCurrentRegion(self, x0, x1, y0, y1, color = 'grey'):
        from matplotlib.patches import Rectangle

        self.ax.add_patch(Rectangle((x0, y0), (x1 - x0), (y1 - y0), alpha=0.1, facecolor = color))
        self.refreshFigure()

    def deleteMarkedPicks(self):
        for shot in self.getShotDict().values():
            for shotnumber in self.getShotsForDeletion():
                if shot.getShotnumber() == shotnumber:
                    for traceID in self.getShotsForDeletion()[shotnumber]:
                        shot.removePick(traceID)
                        print "Deleted the pick for traceID %s on shot number %s" %(traceID, shotnumber)
        self.shots_for_deletion = {} # clear dictionary

    def highlightPicksForShot(self, shot, annotations = False):
        for traceID in shot.getTraceIDlist():
            if shot.getPick(traceID) is not None:
                self.highlightPick(shot, traceID, annotations)
        self.refreshFigure()

    def refreshFigure(self):
        plt.draw()
