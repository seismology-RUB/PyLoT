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
    
    regions.plotTracesInActiveRegions():
     - creates plots (shot.plot_traces) for all traces in the active regions (i.e. chosen by e.g. chooseRectangles)
    
    regions.setAllActiveRegionsForDeletion():
     - highlights all shots in a the active regions for deletion
    
    regions.deleteAllMarkedPicks():
     - deletes the picks (pick flag set to 0) for all shots set for deletion

    regions.deselectSelection(number):
     - deselects the region of number = number
    
    '''

    def __init__(self, ax, cbar, survey):
        self.ax = ax
        self.cbar = cbar
        self.cbv = 'log10SNR'
        self._xlim0 = self.ax.get_xlim()
        self._ylim0 = self.ax.get_ylim()
        self._xlim = self.ax.get_xlim()
        self._ylim = self.ax.get_ylim()
        self.survey = survey
        self.shot_dict = self.survey.getShotDict()
        self._x0 = []
        self._y0 = []
        self._x1 = []
        self._y1 = []
        self._polyx = []
        self._polyy = []
        self.buttons = {}
        self._allpicks = None
        self.shots_found = {}
        self.shots_for_deletion = {}
        self._generateList()
        self._addButtons()
        self.addTextfield()
        self.drawFigure()

    def _generateList(self):
        allpicks = []
        for shot in self.shot_dict.values():
            for traceID in shot.getTraceIDlist():
                allpicks.append((shot.getDistance(traceID),
                 shot.getPickIncludeRemoved(traceID),
                 shot.getShotnumber(),
                 traceID,
                 shot.getFlag(traceID)))

        allpicks.sort()
        self._allpicks = allpicks

    def getShotDict(self):
        return self.shot_dict

    def getShotsForDeletion(self):
        return self.shots_for_deletion

    def _onselect_clicks(self, eclick, erelease):
        '''eclick and erelease are matplotlib events at press and release'''
        print 'region selected x0, y0 = (%3s, %3s), x1, y1 = (%3s, %3s)' % (eclick.xdata,
         eclick.ydata,
         erelease.xdata,
         erelease.ydata)
        x0 = min(eclick.xdata, erelease.xdata)
        x1 = max(eclick.xdata, erelease.xdata)
        y0 = min(eclick.ydata, erelease.ydata)
        y1 = max(eclick.ydata, erelease.ydata)

        shots, numtraces = self.findTracesInShotDict((x0, x1), (y0, y1))
        self.printOutput('Found %d traces in rectangle: %s' % (numtraces, shots))
        key = self.getKey()
        self.shots_found[key] = {'shots': shots,
                                 'selection': 'rect',
                                 'xvalues': (x0, x1),
                                 'yvalues': (y0, y1)}
        self.markRectangle((x0, x1), (y0, y1), key)
        self.disconnectRect()

    def _onselect_verts(self, verts):
        x = verts[0][0]
        y = verts[0][1]
        self._polyx.append(x)
        self._polyy.append(y)

        self.drawPolyLine()

    def _onpress(self, event):
        if event.button == 3:
            self.disconnectPoly()
            self.printOutput('Disconnected polygon selection')

    def addTextfield(self, xpos = 0, ypos = 0.95, width = 1, height = 0.03):
        self.axtext = self.ax.figure.add_axes([xpos,
         ypos,
         width,
         height])
        self.axtext.xaxis.set_visible(False)
        self.axtext.yaxis.set_visible(False)

    def writeInTextfield(self, text = None):
        self.setXYlim(self.ax.get_xlim(), self.ax.get_ylim())
        self.axtext.clear()
        self.axtext.text(0.01, 0.5, text, verticalalignment='center', horizontalalignment='left')
        self.drawFigure()

    def _addButtons(self):
        xpos1 = 0.13
        xpos2 = 0.6
        dx = 0.06
        self.addButton('Rect', self.chooseRectangles, xpos=xpos1, color='white')
        self.addButton('Poly', self.choosePolygon, xpos=xpos1 + dx, color='white')
        self.addButton('Plot', self.plotTracesInActiveRegions, xpos=xpos1 + 2 * dx, color='yellow')
        self.addButton('SNR', self.refreshLog10SNR, xpos=xpos1 + 3 * dx, color='cyan')
        self.addButton('PE', self.refreshPickerror, xpos=xpos1 + 4 * dx, color='cyan')
        self.addButton('SPE', self.refreshSPE, xpos=xpos1 + 5 * dx, color='cyan')
        self.addButton('DesLst', self.deselectLastSelection, xpos=xpos2 + dx, color='green')
        self.addButton('SelAll', self.setAllActiveRegionsForDeletion, xpos=xpos2 + 2 * dx)
        self.addButton('DelAll', self.deleteAllMarkedPicks, xpos=xpos2 + 3 * dx, color='red')

    def addButton(self, name, action, xpos, ypos = 0.91, color = None):
        from matplotlib.widgets import Button
        self.buttons[name] = {'ax': None,
         'button': None,
         'action': action,
         'xpos': xpos}
        ax = self.ax.figure.add_axes([xpos,
         ypos,
         0.05,
         0.03])
        button = Button(ax, name, color=color, hovercolor='grey')
        button.on_clicked(action)
        self.buttons[name]['ax'] = ax
        self.buttons[name]['button'] = button
        self.buttons[name]['xpos'] = xpos

    def getKey(self):
        if self.shots_found.keys() == []:
            key = 1
        else:
            key = max(self.shots_found.keys()) + 1
        return key

    def drawPolyLine(self):
        self.setXYlim(self.ax.get_xlim(), self.ax.get_ylim())
        x = self._polyx
        y = self._polyy
        if len(x) >= 2 and len(y) >= 2:
            self.ax.plot(x[-2:], y[-2:], 'k', alpha=0.1)
        self.drawFigure()

    def drawLastPolyLine(self):
        self.setXYlim(self.ax.get_xlim(), self.ax.get_ylim())
        x = self._polyx
        y = self._polyy
        if len(x) >= 2 and len(y) >= 2:
            self.ax.plot((x[-1], x[0]), (y[-1], y[0]), 'k', alpha=0.1)
        self.drawFigure()

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
        self.printOutput('Found %d traces in polygon: %s' % (numtraces, shots))

    def printOutput(self, text):
        print text
        self.writeInTextfield(text)

    def chooseRectangles(self, event = None):
        '''
        Activates matplotlib widget RectangleSelector.
        '''
        from matplotlib.widgets import RectangleSelector
        if hasattr(self, '_cidPoly'):
            self.disconnectPoly()
        self.printOutput('Select rectangle is active. Press and hold left mousebutton.')
        self._cidRect = None
        self._cidRect = self.ax.figure.canvas.mpl_connect('button_press_event', self._onpress)
        self._rectangle = RectangleSelector(self.ax, self._onselect_clicks)
        return self._rectangle

    def choosePolygon(self, event = None):
        '''
        Activates matplotlib widget LassoSelector.
        '''
        from matplotlib.widgets import LassoSelector
        if hasattr(self, '_cidRect'):
            self.disconnectRect()
        self.printOutput('Select polygon is active. Add points with leftclick. Finish with rightclick.')
        self._cidPoly = None
        self._cidPoly = self.ax.figure.canvas.mpl_connect('button_press_event', self._onpress)
        self._lasso = LassoSelector(self.ax, self._onselect_verts)
        return self._lasso

    def disconnectPoly(self, event = None):
        if not hasattr(self, '_cidPoly'):
            self.printOutput('no poly selection found')
            return
        self.ax.figure.canvas.mpl_disconnect(self._cidPoly)
        del self._cidPoly
        self.finishPolygon()
        self._lasso.disconnect_events()
        print 'disconnected poly selection\n'

    def disconnectRect(self, event = None):
        if not hasattr(self, '_cidRect'):
            self.printOutput('no rectangle selection found')
            return
        self.ax.figure.canvas.mpl_disconnect(self._cidRect)
        del self._cidRect
        self._rectangle.disconnect_events()
        print 'disconnected rectangle selection\n'

    def deselectLastSelection(self, event = None):
        if self.shots_found.keys() == []:
            self.printOutput('No selection found.')
            return
        key = max(self.shots_found.keys())
        self.deselectSelection(key)

    def deselectSelection(self, key, color = 'green', alpha = 0.1):
        if key not in self.shots_found.keys():
            self.printOutput('No selection found.')
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
        self.printOutput('Deselected selection number %d' % key)

    def findTracesInPoly(self, x, y, picks = 'normal', highlight = True):
        def dotproduct(v1, v2):
            return sum((a * b for a, b in zip(v1, v2)))

        def getlength(v):
            return math.sqrt(dotproduct(v, v))

        def getangle(v1, v2):
            return np.rad2deg(math.acos(dotproduct(v1, v2) / (getlength(v1) * getlength(v2))))

        def insidePoly(x, y, pickX, pickY):
            angle = 0
            epsilon = 1e-07
            for index in range(len(x)):
                xval1 = x[index - 1]; yval1 = y[index - 1]
                xval2 = x[index]; yval2 = y[index]
                angle += getangle([xval1 - pickX, yval1 - pickY], [xval2 - pickX, yval2 - pickY])
            if 360 - epsilon <= angle <= 360 + epsilon: ### IMPROVE THAT??
                return True

        if len(x) == 0 or len(y) == 0:
            self.printOutput('No polygon defined.')
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
                        if shotnumber not in shots_found.keys():
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
                if shotnumber not in shots_found.keys():
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

        if shot.getFlag(traceID) is 0:
            return

        self.ax.scatter(shot.getDistance(traceID), shot.getPick(traceID), s = 50, marker = 'o', facecolors = 'none', edgecolors = 'm', alpha = 1)
        if annotations == True:
            self.ax.annotate(s='s%s|t%s' % (shot.getShotnumber(), traceID), xy=(shot.getDistance(traceID), shot.getPick(traceID)), fontsize='xx-small')

    def highlightAllActiveRegions(self):
        '''
        Highlights all picks in all active regions.
        '''
        for key in self.shots_found.keys():
            for shotnumber in self.shots_found[key]['shots'].keys():
                for traceID in self.shots_found[key]['shots'][shotnumber]:
                    self.highlightPick(self.shot_dict[shotnumber], traceID)
        self.drawFigure()

    def plotTracesInActiveRegions(self, event = None, keys = 'all', maxfigures = 20):
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
            self.printOutput('No picks defined in that region(s)')

    def setAllActiveRegionsForDeletion(self, event = None):
        keys = []
        for key in self.shots_found.keys():
            keys.append(key)
        self.setRegionForDeletion(keys)

    def setRegionForDeletion(self, keys):
        if type(keys) == int:
            keys = [keys]

        for key in keys:
            for shotnumber in self.shots_found[key]['shots'].keys():
                if shotnumber not in self.shots_for_deletion:
                    self.shots_for_deletion[shotnumber] = []
                for traceID in self.shots_found[key]['shots'][shotnumber]:
                    if traceID not in self.shots_for_deletion[shotnumber]:
                        self.shots_for_deletion[shotnumber].append(traceID)
            self.deselectSelection(key, color = 'red', alpha = 0.2)

            self.deselectSelection(key, color='red', alpha=0.2)

        self.printOutput('Set region(s) %s for deletion' % keys)

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
        self.ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, alpha=alpha, facecolor=color, linewidth=linewidth))
        if key is not None:
            self.ax.text(x0 + (x1 - x0) / 2, y0 + (y1 - y0) / 2, str(key))
        self.drawFigure()

    def markPolygon(self, x, y, key = None, color = 'grey', alpha = 0.1, linewidth = 1):
        from matplotlib.patches import Polygon
        poly = Polygon(np.array(zip(x, y)), color=color, alpha=alpha, lw=linewidth)
        self.ax.add_patch(poly)
        if key is not None:
            self.ax.text(min(x) + (max(x) - min(x)) / 2, min(y) + (max(y) - min(y)) / 2, str(key))
        self.drawFigure()

    def clearShotsForDeletion(self):
        '''
        Clears the list of shots marked for deletion.
        '''
        self.shots_for_deletion = {}
        print('Cleared all shots that were set for deletion.')

    def getShotsForDeletion(self):
        return self.shots_for_deletion

    def deleteAllMarkedPicks(self, event = None):
        '''
        Deletes all shots set for deletion.
        '''
        if len(self.getShotsForDeletion()) is 0:
            self.printOutput('No shots set for deletion.')
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

    def setXYlim(self, xlim, ylim):
        self._xlim, self._ylim = xlim, ylim

    def refreshLog10SNR(self, event = None):
        cbv = 'log10SNR'
        self.refreshFigure(self, colorByVal=cbv)

    def refreshPickerror(self, event = None):
        cbv = 'pickerror'
        self.refreshFigure(self, colorByVal=cbv)

    def refreshSPE(self, event = None):
        cbv = 'spe'
        self.refreshFigure(self, colorByVal=cbv)

    def refreshFigure(self, event = None, colorByVal = None):
        if colorByVal == None:
            colorByVal = self.cbv
        else:
            self.cbv = colorByVal
        self.printOutput('Refreshing figure...')
        self.ax.clear()
        self.ax = self.survey.plotAllPicks(ax=self.ax, cbar=self.cbar, refreshPlot=True, colorByVal=colorByVal)
        self.setXYlim(self.ax.get_xlim(), self.ax.get_ylim())
        self.markAllActiveRegions()
        self.highlightAllActiveRegions()
        self.drawFigure()
        self.printOutput('Done!')

    def drawFigure(self, resetAxes = True):
        if resetAxes == True:
            self.ax.set_xlim(self._xlim)
            self.ax.set_ylim(self._ylim)
        plt.draw()
