# PyLoT Documentation

- [PyLoT Documentation](#pylot-documentation)
- [PyLoT GUI](#pylot-gui)
  - [First start](#first-start)
  - [Main Screen](#main-screen)
    - [Waveform Plot](#waveform-plot)
      - [Mouse view controls](#mouse-view-controls)
      - [Buttons](#buttons)
    - [Array Map](#array-map)
    - [Eventlist](#eventlist)
  - [Usage](#usage)
    - [Projects and Events](#projects-and-events)
    - [Event folder structure](#event-folder-structure)
    - [Loading event information from CSV file](#loading-event-information-from-csv-file)
    - [Adding events to project](#adding-events-to-project)
    - [Saving projects](#saving-projects)
    - [Adding metadata](#adding-metadata)
- [Picking](#picking)
  - [Manual Picking](#manual-picking)
    - [Picking window](#picking-window)
      - [Picking Window Settings](#picking-window-settings)
    - [Filtering](#filtering)
    - [Export and Import of manual picks](#export-and-import-of-manual-picks)
      - [Export](#export)
      - [Import](#import)
  - [Automatic Picking](#automatic-picking)
    - [Tuning](#tuning)
    - [Production run of the autopicker](#production-run-of-the-autopicker)
    - [Evaluation of automatic picks](#evaluation-of-automatic-picks)
      - [1. Jackknife check](#1-jackknife-check)
      - [2. Wadati check](#2-wadati-check)
    - [Comparison between automatic and manual picks](#comparison-between-automatic-and-manual-picks)
    - [Export and Import of automatic picks](#export-and-import-of-automatic-picks)
- [Location determination](#location-determination)
- [FAQ](#faq)

# PyLoT GUI

This section describes how to use PyLoT graphically to view waveforms and create manual or automatic picks.

## First start

After opening PyLoT for the first time, the setup routine asks for the following information:

Questions:
1. Full Name
2. Authority: Enter authority/institution name
3. Format: Enter output format (*.xml, *.cnv, *.obs)

 [//]: <> (TODO: explain what these things mean, where they are used)

## Main Screen

After entering the [information](#first-start), PyLoTs main window is shown. It defaults to a view of the [Waveform Plot](#waveform-plot), which starts empty.

<img src=images/gui/pylot-main-screen.png alt="Tune autopicks button" title="Tune autopicks button">

Add trace data by [loading a project](#projects-and-events) or by [adding event data](#adding-events-to-project).

### Waveform Plot

The waveform plot shows a trace list of all stations of an event.   
Click on any trace to open the stations [picking window](#picking-window), where you can review automatic and manual picks.

<img src=images/gui/pylot-waveform-plot.png alt="A Waveform Plot showing traces of one event">

Above the traces the currently displayed event can be selected.
In the bottom bar information about the trace under the mouse cursor is shown. This information includes the station name (station), the absolute UTC time (T) of the point under the mouse cursor and the relative time since the first trace start in seconds (t) as well as a trace count.

#### Mouse view controls 

Hold left mouse button and drag to pan view.

Hold right mouse button and
Direction | Result
--- | ---
Move the mouse up | Increase amplitude scale
Move the mouse down | Decrease amplitude scale
Move the mouse right | Increase time scale
Move the mouse left | Decrease time scale

Press right mouse button and click "View All" from the context menu to reset the view.

#### Buttons

[//]: <> (Hack: We need these invisible spaces to add space to the first column, otherwise )

Icon &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; | Description
--- | ---
<img src="../icons/newfile.png" alt="Create new project" width="64" height="64"> | Create a new project, for more information about projects see [Projects and Events](#projects-and-events).
<img src="../icons/openproject.png" alt="Open project" width="64" height="64"> | Load a project file from disk.
<img src="../icons/saveproject.png" alt="Save Project" width="64" height="64"> | Save all current events into an associated project file on disk. If there is no project file currently associated, you will be asked to create a new one.
<img src="../icons/saveprojectas.png" alt="Save Project as" width="64" height="64"> | Save all current events into a new project file on disk. See [Saving projects](#saving-projects).
<img src="../icons/add.png" alt="Add event data" width="64" height="64"> | Add event data by selecting directories containing waveforms. For more information see [Event folder structure](#event-folder-structure).
<img src="../icons/openpick.png" alt="Load event information" width="64" height="64"> | Load picks/origins from disk into the currently displayed event. If a pick already exists for a station, the one from file will overwrite the existing one.
<img src="../icons/openpicks.png" alt="Load information for all events" width="64" height="64"> | Load picks/origins for all events of the current project. PyLoT searches for files within the directory of the event and tries to load them for that event. For this function to work, the files containing picks/origins have to be named as described in [Event folder   structure](#event-folder-structure). If a pick already exists for a station, the one from file will overwrite the existing one.
<img src="../icons/savepicks.png" alt="Save picks" width="64" height="64"> | Save event information such as picks and origin to file. You will be asked to select a directory in which this information should be saved.
<img src="../icons/openloc.png" alt="Load location information" width="64" height="64"> | Load location information from disk,
<img src="../icons/Matlab_PILOT_icon.png" alt="Load legacy information" width="64" height="64"> | Load event information from a previous, MatLab based PILOT version.
<img src="../icons/key_Z.png" alt="Display Z" width="64" height="64"> | Display Z component of streams in waveform plot.
<img src="../icons/key_N.png" alt="Display N" width="64" height="64"> | Display N component of streams in waveform plot.
<img src="../icons/key_E.png" alt="Display E" width="64" height="64"> | Display E component of streams in waveform plot.
<img src="../icons/tune.png" alt="Tune Autopicker" width="64" height="64"> | Open the [Tune Autopicker window](#tuning).
<img src="../icons/autopylot_button.png" alt="" width="64" height="64"> | Opens a window that allows starting the autopicker for all events ([Production run of the AutoPicker](#production-run-of-the-autopicker)).
<img src="../icons/compare_button.png" alt="Comparison" width="64" height="64"> | Compare automatic and manual picks, only available if automatic and manual picks for an event exist. See [Comparison between automatic and manual picks](#comparison-between-automatic-and-manual-picks).
<img src="../icons/locate_button.png" alt="Locate event" width="64" height="64"> | Run a location routine (NonLinLoc) as configured in the settings on the picks. See [Location determination](#location-determination).


### Array Map

The array map will display a color diagram to allow a visual check of the consistency of picks across multiple stations. This works by calculating the time difference of every onset to the earliest onset. Then isolines are drawn between stations with the same time difference and the areas between isolines are colored.  
The result should resemble a color gradient as the wavefront rolls over the network area. Stations where picks are earlier/later than their neighbours can be reviewed by clicking on them, which opens the [picking window](#picking-window).

Above the Array Map the picks that are used to create the map can be customized.
The phase of picks that should be used can be selected, which allows checking the consistency of the P- and S-phase separately.
Additionally the pick type can be set to manual, automatic or hybrid, meaning display only manual picks, automatic picks or only display automatic picks for stations where there are no manual ones.

![Array Map](images/gui/arraymap-example.png "Array Map")  
*Array Map for an event at the Northern Mid Atlantic Ridge, between North Africa and Mexico (Lat. 22.58, Lon. -45.11). The wavefront moved from west to east over the network area (Alps and Balcan region), with the earliest onsets in blue in the west.*

To be able to display an array map PyLoT needs to load an inventory file, where the metadata of seismic stations is kept. For more information see [Metadata](#adding-metadata). Additionally automatic or manual picks need to exist for the current event.

### Eventlist

The eventlist displays event parameters. The displayed parameters are saved in the .xml file in the event folder. Events can be deleted from the project by pressing the red X in the leftmost column of the corresponding event.

<img src="images/gui/eventlist.png" alt="Eventlist">

Column | Description
--- | ---
Event | Full path to the events folder.
Time | Time of event.
Lat | Latitude in degrees of event location.
Lon | Longitude in degrees of event location.
Depth | Depth in km of event.
Mag | Magnitude of event.
[N] MP | Number of manual picks.
[N] AP | Number of automatic picks.
Tuning Set | Select whether this event is a Tuning event. See [Automatic Picking](#automatic-picking).
Test Set | Select whether this event is a Test event. See [Automatic Picking](#automatic-picking).
Notes | Free form text field for notes regarding this event. Text will be saved in the notes.txt file in the event folder.

## Usage

### Projects and Events

PyLoT uses projects to categorize different seismic data. A project consists of one or multiple events. Events contain seismic traces from one or multiple stations. An event also contains further information, e.g. origin time, source parameters and automatic as well as manual picks.
Projects are used to group events which should be analysed together. A project could contain all events from a specific region within a timeframe of interest or all recorded events of a seismological experiment.

### Event folder structure

PyLoT expects the following folder structure for seismic data:
* Every event should be in it's own folder with the following naming scheme for the folders:
   ``e[id].[doy].[yy]``, where ``[id]`` is a four-digit numerical id increasing from 0001, ``[doy]`` the three digit day of year and ``[yy]`` the last two digits of the year of the event. This structure has to be created by the user of PyLoT manually.
* These folders should contain the seismic data for their event as ``.mseed`` or other supported filetype
* All automatic and manual picks should be in an ``.xml`` file in their event folder. PyLoT saves picks in this file. This file does not have to be added manually unless there are picks to be imported. The format used to save picks is QUAKEML.   
Picks are saved in a file with the same filename as the event folder with ``PyLoT_`` prepended.
* The file ``notes.txt`` is used for saving analysts comments. Everything saved here will be displayed in the 'Notes' column of the eventlist.

### Loading event information from CSV file

Event information can be saved in a ``.csv`` file located in the rootpath. The file is made from one header line, which is followed by one or multiple data lines. Values are separated by comma, while a dot is used as a decimal separator.   
This information is then shown in the table in the [Eventlist tab](#Eventlist).

One example header and data line is shown below.
```event,Date,Time,Magnitude,Lat,Long,Depth,Region,Basis Lat,Basis Long,Distance [km],Distance [rad],Distance [deg]```   
```e0001.024.16,24/01/16,10:30:30,7.1,59.66,-153.45,128,Southern Alaska,46.62,10.26,8104.65,1.27,72.89,7.1```

The meaning of the header entries is:

Header | description
--- | ---
event | Event id, has to be the same as the folder name in which waveform data for this event is kept.
Data | Origin date of the event, format DD/MM/YY or DD/MM/YYYY.
Time | Origin time of the event. Format HH:MM:SS.
Lat, Long | Origin latitude and longitude in decimal degrees.
Region | Flinn-Engdahl region name.
Basis Lat, Basis Lon | Latitude and longitude of the basis of the station network in decimal degrees.
Distance [km] | Distance from origin coordinates to basis coordinates in km.
Distance [rad] | Distance from origin coordinates to basis coordinates in rad.


### Adding events to project

PyLoT GUI starts with an empty project. To add events, use the add event data button. Select one or multiple folders containing events. 

[//]: <> (TODO: explain _Directories: Root path, Data path, Database path_)

### Saving projects

Save the current project from the menu with File->Save project or File->Save project as.
PyLoT uses ``.plp`` files to save project information. This file format is not interchangeable between different versions of Python interpreters.
Saved projects contain the automatic and manual picks. Seismic trace data is not included into the ``.plp`` file, but read from its location used when saving the file.

### Adding metadata

[//]: <> (TODO: Add picture of metadata "manager" when it is done)

PyLoT can handle ``.dless``, ``.xml``, ``.resp`` and ``.dseed`` file formats for Metadata. Metadata files stored on disk can be added to a project by clicking *Edit*->*Manage Inventories*. This opens up a window where the folders which contain metadata files can be selected. PyLoT will then search these files for the station names when it needs the information.

# Picking

PyLoTs automatic and manual pick determination works as following:
* Using certain parameters, a first initial/coarse pick is determined. The first manual pick is determined by visual review of the whole waveform and selection of the most likely onset by the analyst. The first automatic pick is determined by calculation of a characteristic function (CF) for the seismic trace. When a wave arrives, the CFs properties change, which is determined as the signals onset.
* Afterwards, a refined set of parameters is applied to a small part of the waveform around the initial onset. For manual picks this means a closer view of the trace, for automatic picks this is done by a recalculated CF with different parameters.
* This second picking phase results in the precise pick, which is treated as the onset time.

## Manual Picking

To create manual picks, you will need to open or create a project that contains seismic trace data (see [Adding events to projects](#adding-events-to-project)). Click on a trace to open the [Picking window](#picking-window).

### Picking window

Open the picking window of a station by leftclicking on any trace in the waveform plot. Here you can create manual picks for the selected station. 

<img src="images/gui/picking/pickwindow.png" alt="Picking window">

*Picking window of a station.*

#### Picking Window Settings

Icon | Shortcut | Menu Alternative | Description
---|---|---|---
<img src="../icons/filter_p.png" alt="Filter P" width="64" height="64"> | p | Filter->Apply P Filter | Filter all channels according to the options specified in Filter parameter, P Filter section.
<img src="../icons/filter_s.png" alt="Filter S" width="64" height="64"> | s | Filter->Apply S Filter | Filter all channels according to the options specified in Filter parameter, S Filter section.
<img src="../icons/key_A.png" alt="Filter Automatically" width="64" height="64"> | Ctrl + a | Filter->Automatic Filtering | If enabled, automatically select the correct filter option (P, S) depending on the selected phase to be picked.
![desc](images/gui/picking/phase_selection.png "Phase selection") | 1 (P) or 5 (S) | Picks->P or S | Select phase to pick. If Automatic Filtering is enabled, this will apply the appropriate filter depending on the phase.
![Zoom into](../icons/zoom_in.png "Zoom into waveform") | - | - | Zoom into waveform.
![Reset zoom](../icons/zoom_0.png "Reset zoom") | - | - | Reset zoom to default view.
![Delete picks](../icons/delete.png "Delete picks") | - | - | Delete all manual picks on this station.
![Rename a phase](../icons/sync.png "Rename a phase") | - | - | Click this button and then the picked phase to rename it.
![Continue](images/gui/picking/continue.png "Continue with next station") | - | - | If checked, after accepting the manual picks for this station with 'OK', the picking window for the next station will be opened. This option is useful for fast manual picking of a complete event.
Estimated onsets | - | - | Show the theoretical onsets for this station. Needs metadata and origin information.
Compare to channel | - | - | Select a data channel to compare against. The selected channel will be displayed in the picking window behind every channel allowing the analyst to visually compare signal correlation between different channels.
Scaling | - | - | Individual means every channel is scaled to its own maximum. If a channel is selected here, all channels will be scaled relatively to this channel.

Menu Command | Shortcut | Description
---|---|---
P Channels and S Channels | - | Select which channels should be treated as P or S channels during picking. When picking a phase, only the corresponding channels will be shown during the precise pick. Normally, the Z channel should be selected for the P phase and the N and E channel for the S phase. 

### Filtering

Access the Filter options by pressing Ctrl+f on the Waveform plot or by the menu under *Edit*->*Filter Parameter*.

<img src=images/gui/pylot-filter-options.png>

 Here you are able to select filter type, order and frequencies for the P and S pick separately. These settings are used in the GUI for displaying the filtered waveform data and during manual picking. The values used by PyLoT for automatic picking are displayed next to the manual values. They can be changed in the [Tune Autopicker dialog](#tuning).  
A green value automatic value means the automatic and manual filter parameter is configured the same, red means they are configured differently.
By toggling the "Overwrite filteroptions" checkmark you can set whether the manual precise/second pick uses the filter settings for the automatic picker (unchecked) or whether it uses the filter options in this dialog (checked).
To guarantee consistent picking results between automatic and manual picking it is recommended to use the same filter settings for the determination of automatic and manual picks.

### Export and Import of manual picks

#### Export

After the creation of manual picks they can either be saved in the project file (see [Saving projects](#saving-projects)). Alternatively the picks can be exported by pressing the <img src="../icons/savepicks.png" alt="Save event information button" title="Save picks button" height=24 width=24> button above the waveform plot or in the menu File->Save event information (shortcut Ctrl+p). Select the event directory in which to save the file. The filename will be ``PyLoT_[event_folder_name].[filetype selected during first startup]``.     
You can rename and copy this file, but PyLoT will then no longer be able to automatically recognize the correct picks for an event and the file will have to be manually selected when loading. 

#### Import

To import previously saved picks press the <img src="../icons/openpick.png" alt="Load event information button" width="24" height="24"> button and select the file to load. You will be asked to save the current state of your project if you have not done so before. You can continue without saving by pressing "Discard". This does not delete any information from your project, it just means that no project file is saved before the changes of importing picks are applied.
PyLoT will automatically load files named after the scheme it uses when saving picks, described in the paragraph above. If it can't find any matching files, a file dialogue will open and you can select the file you wish to load. 

If you see a warning "Mismatch in event identifiers" and are asked whether to continue loading the picks, this means that PyLoT doesn't recognize the picks in the file as belonging to this specific event. They could have either been saved under a different installation of PyLoT but with the same waveform data, which means they are still compatible and you can continue loading them. Or they could be picks from a different event, in which case loading them is not recommended. 

## Automatic Picking


The general workflow for automatic picking is as following:
- After setting up the project by loading waveforms and optionally metadata, the right parameters for the autopicker have to be determined
- This [tuning](#tuning) is done for single stations with immediate graphical feedback of all picking results
- Afterwards the autopicker can be run for all or a subset of events from the project

For automatic picking PyLoT discerns between tune and test events, which the user has to set as such. Tune events are used to calibrate the autopicking algorithm, test events are then used to test the calibration. The purpose of that is controlling whether the parameters found during tuning are able to reliably pick the "unknown" test events.  
If this behaviour is not desired and all events should be handled the same, dont mark any events. Since this is just a way to group events to compare the picking results, nothing else will change.

### Tuning

Tuning describes the process of adjusting the autopicker settings to the characteristics of your data set. To do this in PyLoT, use the <img src=../icons/tune.png height=24 alt="Tune autopicks button" title="Tune autopicks button"> button to open the Tune Autopicker. 

<img src=images/gui/tuning/tune_autopicker.png>

View of a station in the Tune Autopicker window. 
1. Select the event to be displayed and processed.
2. Select the station from the event. 
3. To pick the currently displayed trace, click the <img src=images/gui/tuning/autopick_trace_button.png alt="Pick trace button" title="Autopick trace button" height=16> button.
4. These tabs are used to select the current view. __Traces Plot__ contains a plot of the stations traces, where manual picks can be created/edited. __Overview__ contains graphical results of the automatic picking process. The __P and S tabs__ contain the automatic picking results of the P and S phase, while __log__ contains a useful text output of automatic picking.
5. These buttons are used to load/save/reset settings for automatic picking. The parameters can be saved in PyLoT input files, which have the file ending *.in*. They are human readable text files, which can also be edited by hand. Saving the parameters allows you to load them again later, even on different machines.
6. These menus control the behaviour of the creation of manual picks from the Tune Autopicker window. Picks allows to select the phase for which a manual pick should be created, Filter allows to filter waveforms and edit the filter parameters. P-Channels and S-Channels allow to select the channels that should be displayed when creating a manual P or S pick. 
7. This menu is the same as in the [Picking Window](#picking-window-settings), with the exception of the __Manual Onsets__ options. The __Manual Onsets__ buttons accepts or reject the manual picks created in the Tune Autopicker window, pressing accept adds them to the manual picks for the event, while reject removes them.   
8. The traces plot in the centre allows creating manual picks and viewing the waveforms.
9. The parameters which influence the autopicking result are in the Main settings and Advanced settings tabs on the left side. For a description of all the parameters see the [tuning documentation](tuning.md).

### Production run of the autopicker

After the settings used during tuning give the desired results, the autopicker can be used on the complete dataset. To invoke the autopicker on the whole set of events, click the <img src=../icons/autopylot_button.png alt="Autopick" title="Autopick" height=32> button.

### Evaluation of automatic picks

PyLoT has two internal consistency checks for automatic picks that were determined for an event:
1. Jackknife check
2. Wadati check

#### 1. Jackknife check

The jackknife test in PyLoT checks the consistency of automatically determined P-picks by checking the statistical variance of the picks. The variance of all P-picks is calculated and compared to the variance of subsets, in which one pick is removed.   
The idea is, that picks that are close together in time should not influence the estimation of the variance much, while picks whose positions deviates from the norm influence the variance to a greater extent. If the estimated variance of a subset with a pick removed differs to much from the estimated variance of all picks, the pick that was removed from the subset will be marked as invalid.   
The factor by which picks are allowed to skew from the estimation of variance can be configured, it is called *jackfactor*, see [here](tuning.md#Pick-quality-control).

Additionally, the deviation of picks from the median is checked. For that, the median of all P-picks that passed the Jackknife test is calculated. Picks whose onset times deviate from the mean onset time by more than the *mdttolerance* are marked as invalid.

<img src=images/gui/jackknife_plot.png title="Jackknife/Median test diagram">

*The result of both tests (Jackknife and Median) is shown in a diagram afterwards. The onset time is plotted against a running number of stations. Picks that failed either the Jackknife or the median test are colored red. The median is plotted as a green line.* 

The Jackknife and median check are suitable to check for picks that are outside of the expected time window, for example, when a wrong phase was picked. It won't recognize picks that are in close proximity to the right onset which are just slightly to late/early.

#### 2. Wadati check

The Wadati check checks the consistency of S picks. For this the SP-time, the time difference between S and P onset is plotted against the P onset time. A line is fitted to the points, which minimizes the error. Then the deviation of single picks to this line is checked. If the deviation in seconds is above the *wdttolerance* parameter ([see here](tuning.md#Pick-quality-control)), the pick is marked as invalid.

<img src=images/gui/wadati_plot.png title="Output diagram of Wadati check">

*The Wadati plot in PyLoT shows the SP onset time difference over the P onset time. A first line is fitted (black). All picks which deviate to much from this line are marked invalid (red). Then a second line is fitted which excludes the invalid picks. From this lines slope, the ratio of P and S wave velocity is determined.*

### Comparison between automatic and manual picks

Every pick in PyLoT consists of an earliest possible, latest possible and most likely onset time. 
The earliest and latest possible onset time characterize the uncertainty of a pick. 
This approach is described in Diel, Kissling and Bormann (2012) - Tutorial for consistent phase picking at local to regional distances. 
These times are represented as a Probability Density Function (PDF) for every pick. 
The PDF is implemented as two exponential distributions around the most likely onset as the expected value. 

To compare two single picks, their PDFs are cross correlated to create a new PDF. 
This corresponds to the subtraction of the automatic pick from the manual pick.

 <img src=images/gui/comparison/comparison_pdf.png title="Comparison between automatic and manual pick">

 *Comparison between an automatic and a manual pick for a station in PyLoT by comparing their PDFs.*  
 *The upper plot shows the difference between the two single picks that are shown in the lower plot.*
 *The difference is implemented as a cross correlation between the two PDFs. and results in a new PDF, the comparison PDF.*
 *The expected value of the comparison PDF corresponds to the time distance between the automatic and manual picks most likely onset.*
 *The standard deviation corresponds to the combined uncertainty.*

To compare the automatic and manual picks between multiple stations of an event, the properties of all the comparison PDFs are shown in a histogram.

<img src=images/gui/comparison/compare_widget.png title="Comparison between picks of an event">

*Comparison between the automatic and manual picks for an event in PyLoT.*   
*The top left plot shows the standard deviation of the comparison PDFs for P picks.*
*The bottom left plot shows the expected values of the comparison PDFs for P picks.*
*The top right plot shows the standard deviation of the comparison PDFs for S picks.*
*The bottom right plot shows the expected values of the comparison PDFs for S picks.*   
*The standard deviation plots show that most P picks have an uncertainty between 1 and 2 seconds, while S pick uncertainties have a much larger spread between 1 to 15 seconds.* 
*This means P picks have higher quality classes on average than S picks.*
*The expected values are largely negative, meaning that the algorithm tends to pick earlier than the analyst with the applied settings (Manual - Automatic).*
*The number of samples mentioned in the plots legends is the amount of stations that have an automatic and a manual P pick.*


### Export and Import of automatic picks

Picks can be saved in *.xml* format.

# Location determination

To be added.

# FAQ

Q: During manual picking the error "No channel to plot for phase ..." is displayed, and I am unable to create a pick.  
A: Select a channel that should be used for the corresponding phase in the Pickwindow. For further information read [Picking Window settings](#picking-window-settings).

Q: I see a warning "Mismatch in event identifiers" when loading picks from a file.  
A: This means that PyLoT doesn't recognize the picks in the file as belonging to this specific event. They could have been saved under a different installation of PyLoT but with the same waveform data, which means they are still compatible and you can continue loading them or they could be the picks of a different event, in which case loading them is not recommended.