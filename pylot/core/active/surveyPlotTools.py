# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import math
import numpy as np
plt.interactive(True)

class regions(object):
    '''
    A class used for manual inspection and processing of all picks for the user.

    Examples:

    regions.chooseRectangles():
     - lets the user choose several rectangular regions in the plot

    regions.plotTracesInRegions():
     - creates plots (shot.plot_traces) for all traces in the active regions (i.e. chosen by e.g. chooseRectangles)

    regions.setActiveRegionsForDeletion():
     - highlights all shots in a the active regions for deletion

    regions.deleteMarkedPicks():
     - deletes the picks (pick flag set to 0) for all shots set for deletion

    regions.deselectSelection(number):
     - deselects the region of number = number

    '''
    def __init__(self, ax, survey):
        self.ax = ax
        self.survey = survey
        self.shot_dict = self.survey.getShotDict()
        self._x0 = []
        self._y0 = []
        self._x1 = []
        self._y1 = []
        self._polyx = []
        self._polyy = []
        self._allpicks = None
        self.shots_found = {}
        self.shots_for_deletion = {}
        self._generateList()

    def _onselect_clicks(self, eclick, erelease):
        'eclick and erelease are matplotlib events at press and release'
        #print ' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata)     
        #print ' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata)        
        print 'region selected x0, y0 = (%3s, %3s), x1, y1 = (%3s, %3s)'%(eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata)      
        x0 = min(eclick.xdata, erelease.xdata)
        x1 = max(eclick.xdata, erelease.xdata)
        y0 = min(eclick.ydata, erelease.ydata)
        y1 = max(eclick.ydata, erelease.ydata)

        shots, numtraces = self.findTracesInShotDict((x0, x1), (y0, y1))
        print('Found %d traces in rectangle: %s' %(numtraces, shots))

        key = self.getKey()
        self.shots_found[key] = {'shots': shots,
                                 'selection': 'rect',
                                 'xvalues': (x0, x1),
                                 'yvalues': (y0, y1)}
        self.markRectangle((x0, x1), (y0, y1), key)
        
    def _onselect_verts(self, verts):
        x = verts[0][0]
        y = verts[0][1]
        self._polyx.append(x)
        self._polyy.append(y)

        self.drawPolyLine()

    def _onpress(self, event):
        if event.button == 3:
            self.disconnectPoly()

    def getKey(self):
        if self.shots_found.keys() == []:
            key = 1
        else:
            key = max(self.shots_found.keys()) + 1
        return key

    def drawPolyLine(self):
        x = self._polyx
        y = self._polyy
        if len(x) >= 2 and len(y) >= 2:
            plt.plot(x[-2:], y[-2:], 'k')

    def drawLastPolyLine(self):
        x = self._polyx
        y = self._polyy
        if len(x) >= 2 and len(y) >= 2:
            plt.plot((x[-1], x[0]), (y[-1], y[0]), 'k')

    def finishPolygon(self):
        self.drawLastPolyLine()
        x = self._polyx
        y = self._polyy
        self._polyx = []; self._polyy = []

        key = self.getKey()
        self.markPolygon(x, y, key = key)

        shots, numtraces = self.findTracesInPoly(x, y)
        self.shots_found[key] = {'shots': shots,
                                 'selection': 'poly',
                                 'xvalues': x,
                                 'yvalues': y}

        print('Found %d traces in polygon: %s' %(numtraces, shots))

    def markPolygon(self, x, y, key = None, color = 'grey', alpha = 0.1, linewidth = 1):
        from matplotlib.patches import Polygon
        poly = Polygon(np.array(zip(x, y)), color = color, alpha = alpha, lw = linewidth)
        self.ax.add_patch(poly)
        if key is not None:
            self.ax.text((min(x) + (max(x) - min(x)) / 2), (min(y) + (max(y) - min(y)) / 2), str(key))
        self.drawFigure()

    def disconnectPoly(self):
        self.ax.figure.canvas.mpl_disconnect(self._cid)
        del self._cid
        self.finishPolygon()
        self._lasso.disconnect_events()
        print('disconnected poly selection\n')

    def disconnectRect(self):
        self.ax.figure.canvas.mpl_disconnect(self._cid)
        del self._cid
        self._rectangle.disconnect_events()
        print('disconnected rectangle selection\n')

    def chooseRectangles(self):
        '''
        Activates matplotlib widget RectangleSelector.
        '''
        from matplotlib.widgets import RectangleSelector

        print('Select rectangle is active')
        self._cid = self.ax.figure.canvas.mpl_connect('button_press_event', self._onpress)
        self._rectangle = RectangleSelector(self.ax, self._onselect_clicks)
        return self._rectangle

    def choosePolygon(self):
        '''
        Activates matplotlib widget LassoSelector.
        '''
        from matplotlib.widgets import LassoSelector

        print('Select polygon is active')
        self._cid = self.ax.figure.canvas.mpl_connect('button_press_event', self._onpress)        
        self._lasso = LassoSelector(self.ax, self._onselect_verts)
        return self._lasso
    
    def deselectLastSelection(self):
        if self.shots_found.keys() == []:
            print('No selection found.')
            return
        key = max(self.shots_found.keys())
        self.deselectSelection(key)

    def deselectSelection(self, key, color = 'green', alpha = 0.1):
        if not key in self.shots_found.keys():
            print('No selection found.')
            return
        if color is not None:
            if self.shots_found[key]['selection'] == 'rect':
                self.markRectangle(self.shots_found[key]['xvalues'],
                                   self.shots_found[key]['yvalues'],
                                   key = key, color = color, alpha = alpha,
                                   linewidth = 1)
            elif self.shots_found[key]['selection'] == 'poly':
                self.markPolygon(self.shots_found[key]['xvalues'],
                                 self.shots_found[key]['yvalues'],
                                 key = key, color = color, alpha = alpha,
                                 linewidth = 1)                    
        value = self.shots_found.pop(key)
        print('Deselected selection number %d'% key)
        return

    def _generateList(self):
        allpicks = []
        for shot in self.shot_dict.values():
            for traceID in shot.getTraceIDlist():                           
                allpicks.append((shot.getDistance(traceID), shot.getPickIncludeRemoved(traceID),
                                 shot.getShotnumber(), traceID, shot.getFlag(traceID)))
        allpicks.sort()
        self._allpicks = allpicks

    def getShotDict(self):
        return self.shot_dict

    def getShotsForDeletion(self):
        return self.shots_for_deletion

    def findTracesInPoly(self, x, y, picks = 'normal', highlight = True):
        def dotproduct(v1, v2):
            return sum((a*b) for a, b in zip(v1, v2))

        def getlength(v):
            return math.sqrt(dotproduct(v, v))

        def getangle(v1, v2):
            return np.rad2deg(math.acos(dotproduct(v1, v2) / (getlength(v1) * getlength(v2))))

        def insidePoly(x, y, pickX, pickY):
            angle = 0
            epsilon = 10e-8
            for index in range(len(x)):
                xval1 = x[index - 1]; yval1 = y[index - 1]
                xval2 = x[index]; yval2 = y[index]
                angle += getangle([xval1 - pickX, yval1 - pickY], [xval2 - pickX, yval2 - pickY])
            if 360 - epsilon <= angle <= 360 + epsilon: ### IMPROVE THAT??
                return True

        if len(x) == 0 or len(y) == 0:
            print('No polygon defined.')
            return

        shots_found = {}; numtraces = 0
        x0 = min(x); x1 = max(x)
        y0 = min(y); y1 = max(y)

        shots, numtracesrect = self.findTracesInShotDict((x0, x1), (y0, y1), highlight = False)
        for shotnumber in shots.keys():
            shot = self.shot_dict[shotnumber]
            for traceID in shots[shotnumber]:
                if shot.getFlag(traceID) is not 0:
                    pickX = shot.getDistance(traceID)
                    pickY = shot.getPick(traceID)
                    if insidePoly(x, y, pickX, pickY):
                        if not shotnumber in shots_found.keys():
                            shots_found[shotnumber] = []
                        shots_found[shotnumber].append(traceID)
                        if highlight == True:
                            self.highlightPick(shot, traceID)
                        numtraces += 1

        self.drawFigure()
        return shots_found, numtraces
        
    def findTracesInShotDict(self, (x0, x1), (y0, y1), picks = 'normal', highlight = True):
        '''
        Returns traces corresponding to a certain area in the plot with all picks over the distances.
        '''
        shots_found = {}; numtraces = 0
        if picks == 'normal': pickflag = 0
        elif picks == 'includeCutOut': pickflag = None

        for line in self._allpicks:
            dist, pick, shotnumber, traceID, flag = line
            if flag == pickflag: continue ### IMPROVE THAT
            if (x0 <= dist <= x1 and y0 <= pick <= y1):
                if not shotnumber in shots_found.keys():
                    shots_found[shotnumber] = []
                shots_found[shotnumber].append(traceID)
                if highlight == True:
                    self.highlightPick(self.shot_dict[shotnumber], traceID)
                numtraces += 1

        self.drawFigure()
        return shots_found, numtraces

    def highlightPick(self, shot, traceID, annotations = True):
        '''
        Highlights a single pick for a shot(object)/shotnumber and traceID.
        If annotations == True: Displays shotnumber and traceID in the plot.
        '''
        if type(shot) == int:
            shot = self.survey.getShotDict()[shot]

        self.ax.scatter(shot.getDistance(traceID), shot.getPick(traceID), s = 50, marker = 'o', facecolors = 'none', edgecolors = 'm', alpha = 1)
        if annotations == True:
            self.ax.annotate(s = 's%s|t%s'%(shot.getShotnumber(), traceID), xy = (shot.getDistance(traceID), shot.getPick(traceID)), fontsize = 'xx-small')
        self.ax.set_ylim(shot.getCut())

    def highlightAllRegions(self):
        '''
        Highlights all picks in all active regions.
        '''
        for key in self.shots_found.keys():
            for shotnumber in self.shots_found[key]['shots'].keys():
                for traceID in self.shots_found[key]['shots'][shotnumber]:
                    self.highlightPick(self.shot_dict[shotnumber], traceID)
        self.drawFigure()

    def plotTracesInRegions(self, keys = 'all', maxfigures = 20):
        '''
        Plots all traces in the active region or for all specified keys.

        :param: keys
        :type: int or list

        :param: maxfigures, maximum value of figures opened
        :type: int
        '''        
        count = 0
        if keys == 'all':
            keys = self.shots_found.keys()
        elif type(keys) == int:
            keys = [keys]

        if len(self.shots_found) > 0:
            for shot in self.shot_dict.values():
                for key in keys:
                    for shotnumber in self.shots_found[key]['shots']:
                        if shot.getShotnumber() == shotnumber:
                            for traceID in self.shots_found[key]['shots'][shotnumber]:
                                count += 1
                                if count > maxfigures:
                                    print 'Maximum number of figures (%s) reached. %sth figure was not opened.' %(maxfigures, count)
                                    break
                                shot.plot_traces(traceID)
        else:
            print('No picks defined in that region(s)')

    def setActiveRegionsForDeletion(self):
        keys = []
        for key in self.shots_found.keys():
            keys.append(key)
        self.setRegionForDeletion(keys)

    def setRegionForDeletion(self, keys):
        if type(keys) == int:
            keys = [keys]

        for key in keys:
            for shotnumber in self.shots_found[key]['shots'].keys():
                if not shotnumber in self.shots_for_deletion:
                    self.shots_for_deletion[shotnumber] = []
                for traceID in self.shots_found[key]['shots'][shotnumber]:
                    if not traceID in self.shots_for_deletion[shotnumber]:
                        self.shots_for_deletion[shotnumber].append(traceID)
            self.deselectSelection(key, color = 'red', alpha = 0.2)

        print 'Set region(s) %s for deletion'%keys

    def markAllActiveRegions(self):
        for key in self.shots_found.keys():
            if self.shots_found[key]['selection'] == 'rect':
                self.markRectangle(self.shots_found[key]['xvalues'],
                                   self.shots_found[key]['yvalues'], key = key)
            if self.shots_found[key]['selection'] == 'poly':
                self.markPolygon(self.shots_found[key]['xvalues'],
                                 self.shots_found[key]['yvalues'], key = key)
            

    def markRectangle(self, (x0, x1), (y0, y1), key = None, color = 'grey', alpha = 0.1, linewidth = 1):
        '''
        Mark a rectangular region on the axes.
        '''
        from matplotlib.patches import Rectangle

        self.ax.add_patch(Rectangle((x0, y0), (x1 - x0), (y1 - y0),
                                    alpha = alpha, facecolor = color, linewidth = linewidth))
        if key is not None: 
            self.ax.text((x0 + (x1 - x0) / 2), (y0 + (y1 - y0) / 2), str(key))
        self.drawFigure()

    def refreshFigure(self):
        print('Refreshing figure...')
        self.ax.clear()
        self.ax = self.survey.plotAllPicks(ax = self.ax, refreshPlot = True)
        self.markAllActiveRegions()
        self.drawFigure()
        print('Done!')

    def clearShotsForDeletion(self):
        '''
        Clears the list of shots marked for deletion.
        '''
        self.shots_for_deletion = {}
        print('Cleared all shots that were set for deletion.')

    def getShotsForDeletion(self):
        return self.shots_for_deletion
            
    def deleteMarkedPicks(self):
        '''
        Deletes all shots set for deletion.
        '''
        if len(self.getShotsForDeletion()) is 0:
            print('No shots set for deletion.')
            return

        for shot in self.getShotDict().values():
            for shotnumber in self.getShotsForDeletion():
                if shot.getShotnumber() == shotnumber:
                    for traceID in self.getShotsForDeletion()[shotnumber]:
                        shot.removePick(traceID)
                        print "Deleted the pick for traceID %s on shot number %s" %(traceID, shotnumber)
        self.clearShotsForDeletion()
        self.refreshFigure()

    def highlightPicksForShot(self, shot, annotations = False):
        '''
        Highlight all picks for a given shot.
        '''
        if type(shot) is int:
            shot = self.survey.getShotDict()[shotnumber]

        for traceID in shot.getTraceIDlist():
            if shot.getFlag(traceID) is not 0:
                self.highlightPick(shot, traceID, annotations)
        self.drawFigure()

    def drawFigure(self):
        plt.draw()
