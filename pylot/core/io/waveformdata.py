import logging
import os
from dataclasses import dataclass, field
from typing import Union, List

import numpy as np
from obspy import Stream, read
from obspy.io.sac import SacIOError
from obspy.signal.rotate import rotate2zne

from pylot.core.util.utils import full_range, get_stations


@dataclass
class WaveformData:
    wfdata: Stream = field(default_factory=Stream)
    wforiginal: Union[Stream, None] = None
    wf_alt: Stream = field(default_factory=Stream)
    dirty: bool = False

    def set_wf_data(self, fnames: List[str], fnames_alt: List[str] = None, check_rotated=False, metadata=None, tstart=0, tstop=0):
        self.clear_data()
        fnames = self.check_fname_exists(fnames)
        fnames_alt = self.check_fname_exists(fnames_alt)

        if fnames:
            self.append_wf_data(fnames)
            if fnames_alt:
                self.append_wf_data(fnames_alt, alternative=True)
            self.wfdata, _ = self.check_for_gaps_and_merge(self.wfdata)
            self.check_for_nan(self.wfdata)
            if check_rotated and metadata:
                self.wfdata = self.check4rotated(self.wfdata, metadata, verbosity=0)
            self.trim_station_components(self.wfdata, trim_start=True, trim_end=False)
            self.wforiginal = self.wfdata.copy()
            self.dirty = False
            return True
        return False

    def append_wf_data(self, fnames: List[str], alternative: bool = False):
        data_stream = self.wf_alt if alternative else self.wfdata
        warnmsg = ''
        for fname in set(fnames):
            try:
                data_stream += read(fname)
            except TypeError:
                try:
                    data_stream += read(fname, format='GSE2')
                except Exception as e:
                    try:
                        data_stream += read(fname, format='SEGY')
                    except Exception as e:
                        warnmsg += f'{fname}\n{e}\n'
            except SacIOError as se:
                warnmsg += f'{fname}\n{se}\n'

        if warnmsg:
            print(f'WARNING in appendWFData: unable to read waveform data\n{warnmsg}')

    def clear_data(self):
        self.wfdata = Stream()
        self.wforiginal = None
        self.wf_alt = Stream()

    def reset_wf_data(self):
        if self.wforiginal:
            self.wfdata = self.wforiginal.copy()
        else:
            self.wfdata = Stream()
        self.dirty = False

    def check_fname_exists(self, filenames: List[str]) -> List[str]:
        return [fn for fn in filenames if os.path.isfile(fn)]

    def check_for_gaps_and_merge(self, stream):
        """
        check for gaps in Stream and merge if gaps are found
        :param stream: stream of seismic data
        :type stream: `~obspy.core.stream.Stream`
        :return: data stream, gaps returned from obspy get_gaps
        :rtype: `~obspy.core.stream.Stream`
        """
        gaps = stream.get_gaps()
        if gaps:
            merged = ['{}.{}.{}.{}'.format(*gap[:4]) for gap in gaps]
            stream.merge(method=1)
            print('Merged the following stations because of gaps:')
            for merged_station in merged:
                print(merged_station)

        return stream, gaps

    def check_for_nan(self, stream):
        """
        Replace all NaNs in data with nan_value (in place)
        :param stream: stream of seismic data
        :type stream: `~obspy.core.stream.Stream`
        :param nan_value: value which all NaNs are set to
        :type nan_value: float, int
        :return: None
        """
        if not stream:
            return
        for trace in stream:
            np.nan_to_num(trace.data, copy=False, nan=0.)


    def check4rotated(self, stream, metadata=None, verbosity=1):
        """
        Check all traces in data. If a trace is not in ZNE rotation (last symbol of channel code is numeric) and the trace
        is in the metadata with azimuth and dip, rotate it to classical ZNE orientation.
        Rotating the traces requires them to be of the same length, so, all traces will be trimmed to a common length as a
        side effect.
        :param stream: stream object containing seismic traces
        :type stream: `~obspy.core.stream.Stream`
        :param metadata: tuple containing metadata type string and metadata parser object
        :type metadata: (str, `~obspy.io.xseed.parser.Parser`)
        :param verbosity: if 0 print no information at runtime
        :type verbosity: int
        :return: stream object with traditionally oriented traces (ZNE) for stations that had misaligned traces (123) before
        :rtype: `~obspy.core.stream.Stream`
        """

        def rotation_required(trace_ids):
            """
            Derive if any rotation is required from the orientation code of the input.

            :param trace_ids: string identifier of waveform data trace
            :type trace_ids: List(str)
            :return: boolean representing if rotation is necessary for any of the traces
            :rtype: bool
            """
            orientations = [trace_id[-1] for trace_id in trace_ids]
            return any([orientation.isnumeric() for orientation in orientations])

        def rotate_components(wfs_in, metadata=None):
            """
            Rotate components if orientation code is numeric (= non traditional orientation).

            Azimut and dip are fetched from metadata. To be rotated, traces of a station have to be cut to the same length.
            Returns unrotated traces of no metadata is provided
            :param wfs_in: stream containing seismic traces of a station
            :type wfs_in: `~obspy.core.stream.Stream`
            :param metadata: tuple containing metadata type string and metadata parser object
            :type metadata: (str, `~obspy.io.xseed.parser.Parser`)
            :return: stream object with traditionally oriented traces (ZNE)
            :rtype: `~obspy.core.stream.Stream`
            """

            if len(wfs_in) < 3:
                print(f"Stream {wfs_in=}, has not enough components to rotate.")
                return wfs_in

            # check if any traces in this station need to be rotated
            trace_ids = [trace.id for trace in wfs_in]
            if not rotation_required(trace_ids):
                logging.debug(f"Stream does not need any rotation: Traces are {trace_ids=}")
                return wfs_in

            # check metadata quality
            t_start = full_range(wfs_in)
            try:
                azimuths = []
                dips = []
                for tr_id in trace_ids:
                    azimuths.append(metadata.get_coordinates(tr_id, t_start)['azimuth'])
                    dips.append(metadata.get_coordinates(tr_id, t_start)['dip'])
            except (KeyError, TypeError) as err:
                logging.error(
                    f"{type(err)=} occurred: {err=} Rotating not possible, not all azimuth and dip information "
                    f"available in metadata. Stream remains unchanged.")
                return wfs_in
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                raise

            # to rotate all traces must have same length, so trim them
            wfs_out = self.trim_station_components(wfs_in, trim_start=True, trim_end=True)
            try:
                z, n, e = rotate2zne(wfs_out[0], azimuths[0], dips[0],
                                     wfs_out[1], azimuths[1], dips[1],
                                     wfs_out[2], azimuths[2], dips[2])
                print('check4rotated: rotated trace {} to ZNE'.format(trace_ids))
                # replace old data with rotated data, change the channel code to ZNE
                z_index = dips.index(min(
                    dips))  # get z-trace index, z has minimum dip of -90 (dip is measured from 0 to -90, with -90
                # being vertical)
                wfs_out[z_index].data = z
                wfs_out[z_index].stats.channel = wfs_out[z_index].stats.channel[0:-1] + 'Z'
                del trace_ids[z_index]
                for trace_id in trace_ids:
                    coordinates = metadata.get_coordinates(trace_id, t_start)
                    dip, az = coordinates['dip'], coordinates['azimuth']
                    trace = wfs_out.select(id=trace_id)[0]
                    if az > 315 or az <= 45 or 135 < az <= 225:
                        trace.data = n
                        trace.stats.channel = trace.stats.channel[0:-1] + 'N'
                    elif 45 < az <= 135 or 225 < az <= 315:
                        trace.data = e
                        trace.stats.channel = trace.stats.channel[0:-1] + 'E'
            except ValueError as err:
                print(f"{err=} Rotation failed. Stream remains unchanged.")
                return wfs_in

            return wfs_out

        if metadata is None:
            if verbosity:
                msg = 'Warning: could not rotate traces since no metadata was given\nset Inventory file!'
                print(msg)
            return stream
        stations = get_stations(stream)
        for station in stations:  # loop through all stations and rotate data if neccessary
            wf_station = stream.select(station=station)
            rotate_components(wf_station, metadata)
        return stream

    def trim_station_components(stream, trim_start=True, trim_end=True):
        """
        cut a stream so only the part common to all three traces is kept to avoid dealing with offsets
        :param stream: stream of seismic data
        :type stream: `~obspy.core.stream.Stream`
        :param trim_start: trim start of stream
        :type trim_start: bool
        :param trim_end: trim end of stream
        :type trim_end: bool
        :return: data stream
        :rtype: `~obspy.core.stream.Stream`
        """
        starttime = {False: None}
        endtime = {False: None}

        stations = get_stations(stream)

        print('trim_station_components: Will trim stream for trim_start: {} and for '
              'trim_end: {}.'.format(trim_start, trim_end))
        for station in stations:
            wf_station = stream.select(station=station)
            starttime[True] = max([trace.stats.starttime for trace in wf_station])
            endtime[True] = min([trace.stats.endtime for trace in wf_station])
            wf_station.trim(starttime=starttime[trim_start], endtime=endtime[trim_end])

        return stream
