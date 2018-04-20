# Table of contents

- [PyLoT GUI](#pylot-gui)
  * [First start](#first-start)
  * [Main Screen](#main-screen)
    + [Waveform Plot](#waveform-plot)
      - [Mouse controls :](#mouse-controls--)
    + [Array Map](#array-map)
    + [Eventlist](#eventlist)
  * [Usage](#usage)
    + [Projects](#projects)
    + [Event folder structure](#event-folder-structure)
    + [Adding events to project](#adding-events-to-project)
    + [Saving projects](#saving-projects)
- [Manual Picking](#manual-picking)
  * [Picking window](#picking-window)
    + [Picking Window Settings](#picking-window-settings)
  * [Filtering](#filtering)
  * [Export of manual picks](#export-of-manual-picks)
- [Automatic Picking](#automatic-picking)
  * [Tuning](#tuning)
  * [Production run of the autopicker](#production-run-of-the-autopicker)
  * [Evaluation of automatic picks](#evaluation-of-automatic-picks)
  * [Export of automatic picks](#export-of-automatic-picks)

# PyLoT GUI

## First start

Questions:
1. Full Name
2. Authority: Enter authority/instituiton name
3. Format: Enter output format (*.xml, *.cnv, *.obs)

## Main Screen

Add trace data by [loading a project](#projects) or by [adding event data](#adding-events-to-projects).

### Waveform Plot

Click on any trace to open the stations picking window.
In the bottom bar the station name (station), the absolute UTC time (T) of the point under the mouse cursor and the relative time since the first trace start in seconds (t) is shown.

#### Mouse view controls : 

Hold left mouse button and drag to pan view.

Hold right mouse button and
Direction | Result
--- | ---
Move the mouse up | Increase amplitude scale
Move the mouse down | Decrease amplitude scale
Move the mouse right | Increase time scale
Move the mouse left | Decrease time scale

Press right mouse button and click "View All" from the context menu to reset the view.

### Array Map

The array map will display a color diagram to allow checking the consistency of picks across multiple stations.

![Array Map](images/gui/arraymap.png "Array Map")

To be able to display an array map PyLoT needs to load an inventory file, where the metadata of seismic stations is kept. Possible file types are ``.dless``, ``.xml``, ``.resp`` and ``.dseed``.

### Eventlist

The eventlist displays event parameters. The displayed parameters are saved in the .xml file in the event folder. Events can be deleted from the project by pressing the red X in the leftmost column of the corresponding event.

## Usage

### Projects

PyLoT uses projects to categorize different seismic data. A project consists of one or multiple events. Events contain seismic traces from one or multiple stations. An event also contains further information, e.g. origin time, source parameters and automatic and manual picks.
Projects are used to group events which should be analyzed together. A project could contain all events from a specific region within a timeframe of interest or all recorded events of a seismological experiment.

### Event folder structure

PyLoT expects the following folder structure for seismic data:
* Every event should be in it's own folder with the following naming scheme for the folders:
   ``e[id].[doy].[yy]``, where ``[id]`` is a four-digit numerical id increasing from 0001, ``[doy]`` the three digit day of year and ``[yy]`` the last two digits of the year of the event. This structure has to be created by the user of PyLoT manually.
* These folders should contain the seismic data for their event as ``.mseed`` or other supported filetype TODO: LINK HERE
* All automatic and manual picks should be in an ``.xml`` file in their event folder. PyLoT saves picks in this file. This file does not have to be added manually unless there are picks to be imported. The format used to save picks is QUAKEML.
* The file ``notes.txt`` is used for saving analysts comments. Everything saved here will be displayed in the 'Notes' column of the eventlist.

### Adding events to project

PyLoT GUI starts with an empty project. To add events, use the add event data button. Select one or multiple folders containing events. 
TODO: explain _Directories: Root path, Data path, Database path_

### Saving projects

Save the current project from the menu with File->Save project or File->Save project as.
PyLoT uses ``.plp`` files to save project information. This file format is not interchangable between different versions of Python interpreters.
Saved projects contain the automatic and manual picks. Seismic trace data is not included into the ``.plp`` file, but read from its location used when saving the file.

# Manual Picking

To create manual picks, you will need to open or create a project that contains seismic trace data (see [Adding events to projects](#adding-events-to-project)). Click on a trace to open the [Picking window](#picking-window).

## Picking window

Open the picking window of a station by leftclicking on any trace in the waveform plot. Here you can create manual picks for the selected station.

### Picking Window Settings

Icon | Shortcut | Menu Alternative | Description
---|---|---|---
<img src="../icons/filter_p.png" alt="Filter P" width="64" height="64"> | 1 or p | Filter->Apply P Filter | Filter all channels according to the options specified in Filter parameter, P Filter section
<img src="../icons/filter_s.png" alt="Filter S" width="64" height="64"> | 5 or s | Filter->Apply S Filter | Filter all channels according to the options specified in Filter parameter, S Filter section
<img src="../icons/key_A.png" alt="Filter Automatically" width="64" height="64"> | Ctrl + a | Filter->Automatic Filtering | Select the correct filter option (P, S) depending on the selected phase to be picked
<img src="images/gui/picking/phase_selection.png" alt="Phase selection" > | 1 or 5 | Picks->P or S | Select phase to pick
![Zoom into](../icons/zoom_in.png "Zoom into waveform") | - | - | Zoom into waveform TODO: expand
![Reset zoom](../icons/zoom_0.png "Reset zoom") | - | - | Reset zoom to default view
![Delete picks](../icons/delete.png "Delete picks") | - | - | Delete all manual picks on this station
![Rename a phase](../icons/sync.png "Rename a phase") | - | - | Click this button and then the picked phase to rename it
![Continue](images/gui/picking/continue.png "Continue") | - | - | If checked, after accepting the manual picks for this station with 'OK'. the picking window for the next station will be opened
Estimated onsets | - | - | Show the theoretical onsets for this station
Compare to channel | - | - | Select a data channel to compare against. The selected channel will be displayed in the picking window behind every channel to compare signal correlation
Scale by | - | - | Normalized means every channel is scaled to it maximum seperately. If a channel is selected here, all channels will be scaled with regards to this channel
 
Menu Command | Shortcut | Description
---|---|---
P Channels and S Channels | - | Select which channels should be treated as P or S channels during picking. When picking a phase, only the corresponding channels will be shown during the precise pick

## Filtering

Access the Filter options by pressing Ctrl+f on the Waveform plot or by the menu under Edit->Filter Parameter. Here you are able to select filter type, order and frequencies for the P and S pick seperately. These settings are used in the GUI for filtering during manual picking. The values used by PyLoT for automatic picking are displayed next to the manual values.
By toggling the "Overwrite filteroptions" checkmark you can set whether the manual second pick uses the filter settings for the automatic picker (unchecked) or whether to use the filter options in this dialog (checked).
To guarantee consistent picking results between automatic and manual picking it is recommended to use the same fiter settings for the determination of automatic and manual picks.

## Export and Import of manual picks

### Export

After the creation of manual picks they can either be save in the project file (see [Saving projects](#saving-projects)). Alternatively the picks can be exported by pressing the <img src="../icons/savepicks.png" alt="Save event information button" title="Save picks button" height=24 width=24> button above the waveform plot or in the menu File->Save event information (shortcut Ctrl+p). Select the event directory in which to save the file. The filename will be ``PyLoT_[event_folder_name].[filetype selected during first startup]``. TODO: Is the filetype during first startup selected for saving files?

You can rename and copy this file, but PyLoT will then no longer be able to automatically recognize the correct picks for an event and the file will have to be manually selected when loading. 

### Import

To import previously saved picks press the <img src="../icons/openpick.png" alt="Load event information button" width="24" height="24"> button and select the file to load. You are asked to save the current state of your current project if you have not done so before. You can continue without saving by pressing "Discard".
 PyLoT will automaticall load files named after the scheme it uses when saving picks, described in the paragraph above. If it cant file any aptly named files, a file dialog will open and you can select the file you wish to load. 

If you see a warning "Mismatch in event identifiers" and are asked whether to continue loading the picks, this means that PyLoT doesn't recognize the picks in the file as belonging to this specific event. They could have either been saved under a different installation of PyLoT but with the same waveform data, which means they are still compatible and you can continue loading them. Or they could be picks from a different event, in which case loading them is not reccommended. 

# Automatic Picking

## Tuning

To adjust the autopicker settings to the characteristics of your data set, use the <img src=../icons/tune.png height=24 alt="Tune autopicks button" title="Tune autopicks button"> button to open the Tuning Dialog. In the right hand side of the window the Main Settings and Advanced Settings control the result of the automatic picking. To pick the currently displayed trace, click the <img src=images/gui/tuning/autopick_trace_button.png alt="Pick trace button" title="Autopick trace button" height=16> button in the top right corner. 

## Production run of the autopicker

## Evaluation of automatic picks

## Export and Import of automatic picks

# FAQ

Q: During manual picking the error "No channel to plot for phase ..." is displayed, and I am unable to create a pick.  
A: Select a channel that should be used for the corresponding phase in the Pickwindow. For further information read [Picking Window settings](#picking-window-settings).

Q: I see a warning "Mismatch in event identifiers" when loading picks from a file.  
A: This means that PyLoT doesn't recognize the picks in the file as belonging to this specific event. They could have been saved under a different installation of PyLoT but with the same waveform data, which means they are still compatible and you can continue loading them or they could be the picks of a different event, in which case loading them is not reccommended
