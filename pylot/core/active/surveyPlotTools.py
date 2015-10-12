import matplotlib.pyplot as plt
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
        self._allpicks = None
        self.shots_found = {}
        self.shots_for_deletion = {}
        self._generateList()

    def _onselect(self, eclick, erelease):
        'eclick and erelease are matplotlib events at press and release'
        #print ' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata)     
        #print ' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata)        
        print 'region selected x0, y0 = (%3s, %3s), x1, y1 = (%3s, %3s)'%(eclick.xdata, eclick.ydata, erelease.xdata, erelease.ydata)      
        x0 = min(eclick.xdata, erelease.xdata)
        x1 = max(eclick.xdata, erelease.xdata)
        y0 = min(eclick.ydata, erelease.ydata)
        y1 = max(eclick.ydata, erelease.ydata)

        shots = self.findTracesInShotDict((x0, x1), (y0, y1))
        if self.shots_found.keys() == []:
            key = 1
        else:
            key = max(self.shots_found.keys()) + 1

        self.shots_found[key] = {'shots': shots,
                                 'distbin': (x0, x1),
                                 'pickbin': (y0, y1)}
        self.markRegion((x0, x1), (y0, y1), key)
        
    def chooseRectangles(self):
        '''
        Activates matplotlib widget RectangleSelector.
        '''
        from matplotlib.widgets import RectangleSelector

        print 'Select rectangle is active'
        return RectangleSelector(self.ax, self._onselect)

    def deselectLastSelection(self):
        if self.shots_found.keys() == []:
            print('No selection found.')
            return
        key = max(self.shots_found.keys())
        self.deselectSelection(key)

    def deselectSelection(self, key, color = 'green', alpha = 0.1):
        try:
            if color is not None:
                self.markRegion(self.shots_found[key]['distbin'],
                                self.shots_found[key]['pickbin'],
                                key = key, color = color, alpha = alpha, linewidth = 0)
            value = self.shots_found.pop(key)
            print('Deselected selection number %d'% key)
            return
        except:
            print('No selection found.')
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

    def findTracesInShotDict(self, (x0, x1), (y0, y1), picks = 'normal'):
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
                numtraces += 1

        print('Found %d traces: %s' %(numtraces, shots_found))
        return shots_found

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
            self.markRegion(self.shots_found[key]['distbin'],
                            self.shots_found[key]['pickbin'], key = key)
            

    def markRegion(self, (x0, x1), (y0, y1), key = None, color = 'grey', alpha = 0.1, linewidth = 0.1):
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
