import unittest
from unittest import skip
import obspy
from obspy import UTCDateTime
import os
import sys
from pylot.core.pick.autopick import autopickstation
from pylot.core.io.inputs import PylotParameter
from pylot.core.io.data import Data
from pylot.core.util.utils import trim_station_components


class HidePrints:
    """
    Used to hide all standard output the Function to be tested have, since it clutters the test results.
    """

    def __init__(self, hide_prints=True):
        """Create object with hide_prints=False to disable print hiding"""
        self.hide = hide_prints

    def __enter__(self):
        if self.hide:
            self._original_stdout = sys.stdout
            devnull = open(os.devnull, "w")
            sys.stdout = devnull

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.hide:
            sys.stdout = self._original_stdout


class MockMetadata:
    """Mock metadata object used for taupy to avoid reading large dless file from disk.
    get_coordinates must take the same arguments as pylot.core.utils.dataprocssing.py/class Metadata."""

    def __init__(self):
        self.station_names = ['GR.GRA1', 'GR.GRA2', 'G.ECH', 'CH.FIESA', 'Z3.A106A']
        gra1 = {u'azimuth': 0.0, u'dip': -90.0, u'elevation': 499.5, u'latitude': 49.691888, u'local_depth': 0.0,
                u'longitude': 11.22172}
        gra2 = {u'azimuth': 0.0, u'dip': -90.0, u'elevation': 512.0, u'latitude': 49.655208, u'local_depth': 0.0,
                u'longitude': 11.359444}
        ech = {u'azimuth': 90.0, u'dip': 0.0, u'elevation': 580.0, u'latitude': 48.216313, u'local_depth': 250.0,
               u'longitude': 7.158961}
        fiesa = {'azimuth': 0.0, 'dip': -90.0, 'elevation': 2340.5, 'latitude': 46.43521, 'local_depth': 0.0,
                 'longitude': 8.11051}
        a106 = {'azimuth': 90.0, 'dip': 0.0, 'elevation': 468.0, 'latitude': 48.753388, 'local_depth': 0.0,
                'longitude': 9.721937}

        self.coordinates = [gra1, gra2, ech, fiesa, a106]

    def get_coordinates(self, station_id):
        """
        Mocks the method get_coordinates from obspy.io.xseed.parser.Parser object
        to avoid building a parser for the unit tests
        :param station_id: 'GR.GRA1..LHZ' or similar
        :type station_id: str
        :return: dictionary containing azimuth, dip, elevation, latitude, longitude,
        local depth as keys
        :rtype: dict

        >>>m = MockMetadata(); m.get_coordinates('GR.GRA2..LHZ')
        {u'azimuth': 0.0, u'dip': -90.0, u'elevation': 512.0, u'latitude': 49.655208, u'local_depth': 0.0, u'longitude': 11.359444}
        """

        for index, name in enumerate(self.station_names):
            if station_id.startswith(name):
                return self.coordinates[index]


class TestAutopickStation(unittest.TestCase):
    """
    Test the autopickstation function as if it were called from GUI.
    Three stations (GR.GRA1, GR.GRA2, G.ECH) are tested with and without TauPy respectively
    """

    def setUp(self):
        self.event_id = 'e0001.024.16'
        # Create wfstream for picking
        mseed_relative_path = os.path.join(os.path.dirname(__file__), self.event_id, '*.mseed')
        self.wfstream = obspy.read(mseed_relative_path)
        # trim waveform to get the same results as the GUI call
        with HidePrints():
            self.wfstream = trim_station_components(self.wfstream, trim_start=True, trim_end=False)
        self.gra1 = self.wfstream.select(station='GRA1')
        self.gra2 = self.wfstream.select(station='GRA2')
        self.ech = self.wfstream.select(station='ECH')
        self.fiesa = self.wfstream.select(station='FIESA')
        self.a106 = self.wfstream.select(station='A106A')
        # Create input parameter container
        self.inputfile_taupy_enabled = os.path.join(os.path.dirname(__file__), 'autoPyLoT_global_taupy_true.in')
        self.inputfile_taupy_disabled = os.path.join(os.path.dirname(__file__), 'autoPyLoT_global_taupy_false.in')
        self.pickparam_taupy_enabled = PylotParameter(fnin=self.inputfile_taupy_enabled)
        self.pickparam_taupy_disabled = PylotParameter(fnin=self.inputfile_taupy_disabled)
        self.xml_file = os.path.join(os.path.dirname(__file__),self.event_id, 'PyLoT_'+self.event_id+'.xml')
        self.data = Data(evtdata=self.xml_file)
        # create origin for taupy testing
        self.origin = [obspy.core.event.origin.Origin(magnitude=7.1, latitude=59.66, longitude=-153.45, depth=128.0, time=UTCDateTime("2016-01-24T10:30:30.0"))]
        # mocking metadata since reading it takes a long time to read from file
        self.metadata = MockMetadata()

        # show complete diff when difference in results dictionaries are found
        self.maxDiff

    #@skip("Works")
    def test_autopickstation_taupy_disabled_gra1(self):
        expected = {'P': {'picker': 'auto', 'snrdb': 15.405649120980094, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 34.718816470730317, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 31, 690000), 'w0': None, 'spe': 0.93333333333333235, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 28, 890000), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 32, 690000), 'fm': 'D', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.669661906545489, 'network': u'GR', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 30, 690000), 'snr': 11.667187857573905, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 21, 690000), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 29, 690000), 'fm': None, 'spe': 2.6666666666666665, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.gra1, pickparam=self.pickparam_taupy_disabled, metadata=(None, None))
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('GRA1', station)

    def test_autopickstation_taupy_enabled_gra1(self):
        expected = {'P': {'picker': 'auto', 'snrdb': 15.599905299126778, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 36.307013769185403, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 27, 690000), 'w0': None, 'spe': 0.93333333333333235, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 24, 890000), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 28, 690000), 'fm': 'U', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.669661906545489, 'network': u'GR', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 30, 690000), 'snr': 11.667187857573905, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 21, 690000), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 29, 690000), 'fm': None, 'spe': 2.6666666666666665, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.gra1, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin=self.origin)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('GRA1', station)

    def test_autopickstation_taupy_disabled_gra2(self):
        expected = {'P': {'picker': 'auto', 'snrdb': None, 'weight': 9, 'Mo': None, 'marked': 'shortsignallength', 'Mw': None, 'fc': None, 'snr': None, 'mpp': UTCDateTime(2016, 1, 24, 10, 36, 59, 150000), 'w0': None, 'spe': None, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 36, 43, 150000), 'lpp': UTCDateTime(2016, 1, 24, 10, 37, 15, 150000), 'fm': 'N', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': u'GR', 'weight': 4, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 37, 15, 150000), 'snr': None, 'epp': UTCDateTime(2016, 1, 24, 10, 36, 43, 150000), 'mpp': UTCDateTime(2016, 1, 24, 10, 36, 59, 150000), 'fm': None, 'spe': None, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.gra2, pickparam=self.pickparam_taupy_disabled, metadata=(None, None))
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('GRA2', station)

    def test_autopickstation_taupy_enabled_gra2(self):
        expected = {'P': {'picker': 'auto', 'snrdb': 13.957959025719253, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 24.876879503607871, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 29, 150000), 'w0': None, 'spe': 1.0, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 26, 150000), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 30, 150000), 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.573236990555648, 'network': u'GR', 'weight': 1, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 34, 150000), 'snr': 11.410999834108294, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 21, 150000), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 33, 150000), 'fm': None, 'spe': 4.666666666666667, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.gra2, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin = self.origin)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('GRA2', station)

    def test_autopickstation_taupy_disabled_ech(self):
        expected = {'P': {'picker': 'auto', 'snrdb': None, 'weight': 9, 'Mo': None, 'marked': 'SinsteadP', 'Mw': None, 'fc': None, 'snr': None, 'mpp': UTCDateTime(2016, 1, 24, 10, 26, 57), 'w0': None, 'spe': None, 'network': u'G', 'epp': UTCDateTime(2016, 1, 24, 10, 26, 41), 'lpp': UTCDateTime(2016, 1, 24, 10, 27, 13), 'fm': 'N', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': u'G', 'weight': 4, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 27, 13), 'snr': None, 'epp': UTCDateTime(2016, 1, 24, 10, 26, 41), 'mpp': UTCDateTime(2016, 1, 24, 10, 26, 57), 'fm': None, 'spe': None, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.ech, pickparam=self.pickparam_taupy_disabled)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('ECH', station)

    def test_autopickstation_taupy_enabled_ech(self):
        # this station has a long time of before the first onset, so taupy will help during picking
        expected = {'P': {'picker': 'auto', 'snrdb': 9.9753586609166316, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 9.9434218804137107, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 34), 'w0': None, 'spe': 1.6666666666666667, 'network': u'G', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 29), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 35), 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 12.698999454169567, 'network': u'G', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 44), 'snr': 18.616581906366577, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 33), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 43), 'fm': None, 'spe': 3.3333333333333335, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.ech, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin=self.origin)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('ECH', station)

    def test_autopickstation_taupy_disabled_fiesa(self):
        # this station has a long time of before the first onset, so taupy will help during picking
        expected = {'P': {'picker': 'auto', 'snrdb': None, 'weight': 9, 'Mo': None, 'marked': 'SinsteadP', 'Mw': None, 'fc': None, 'snr': None, 'mpp': UTCDateTime(2016, 1, 24, 10, 35, 58), 'w0': None, 'spe': None, 'network': u'CH', 'epp': UTCDateTime(2016, 1, 24, 10, 35, 42), 'lpp': UTCDateTime(2016, 1, 24, 10, 36, 14), 'fm': 'N', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': u'CH', 'weight': 4, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 36, 14), 'snr': None, 'epp': UTCDateTime(2016, 1, 24, 10, 35, 42), 'mpp': UTCDateTime(2016, 1, 24, 10, 35, 58), 'fm': None, 'spe': None, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.fiesa, pickparam=self.pickparam_taupy_disabled)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('FIESA', station)

    def test_autopickstation_taupy_enabled_fiesa(self):
        # this station has a long time of before the first onset, so taupy will help during picking
        expected = {'P': {'picker': 'auto', 'snrdb': 13.921049277904373, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 24.666352170589487, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 47), 'w0': None, 'spe': 1.2222222222222285, 'network': u'CH', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 43, 333333), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 48), 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.893086316477728, 'network': u'CH', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 51, 5), 'snr': 12.283118216397849, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 59, 333333), 'mpp': UTCDateTime(2016, 1, 24, 10, 51, 2), 'fm': None, 'spe': 2.8888888888888764, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.fiesa, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin=self.origin)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual('FIESA', station)

    def test_autopickstation_gra1_z_comp_missing(self):
        """Picking on a stream without a vertical trace should return None"""
        wfstream = self.gra1.copy()
        wfstream = wfstream.select(channel='*E') + wfstream.select(channel='*N')
        with HidePrints():
            result, station = autopickstation(wfstream=wfstream, pickparam=self.pickparam_taupy_disabled, metadata=(None, None))
        self.assertIsNone(result)
        self.assertEqual('GRA1', station)

    def test_autopickstation_gra1_horizontal_comps_missing(self):
        """Picking on a stream without horizontal traces should still pick the P phase on the vertical component"""
        wfstream = self.gra1.copy()
        wfstream = wfstream.select(channel='*Z')
        expected =  {'P': {'picker': 'auto', 'snrdb': 15.405649120980094, 'network': u'GR', 'weight': 0, 'Ao': None, 'Mo': None, 'marked': [], 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 32, 690000), 'Mw': None, 'fc': None, 'snr': 34.718816470730317, 'epp': UTCDateTime(2016, 1, 24, 10, 41, 28, 890000), 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 31, 690000), 'w0': None, 'spe': 0.9333333333333323, 'fm': 'D', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': None, 'weight': 4, 'Mo': None, 'Ao': None, 'lpp': None, 'Mw': None, 'fc': None, 'snr': None, 'marked': [], 'mpp': None, 'w0': None, 'spe': None, 'epp': None, 'fm': 'N', 'channel': None}}
        with HidePrints():
            result, station = autopickstation(wfstream=wfstream, pickparam=self.pickparam_taupy_disabled, metadata=(None, None))
        self.assertEqual(expected, result)
        self.assertEqual('GRA1', station)

    def test_autopickstation_a106_taupy_enabled(self):
        """This station has invalid values recorded on both N and E component, but a pick can still be found on Z"""
        expected = {'P': {'picker': 'auto', 'snrdb': 12.862128789922826, 'network': u'Z3', 'weight': 0, 'Ao': None, 'Mo': None, 'marked': [], 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 34), 'Mw': None, 'fc': None, 'snr': 19.329155459132608, 'epp': UTCDateTime(2016, 1, 24, 10, 41, 30), 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 33), 'w0': None, 'spe': 1.6666666666666667, 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': u'Z3', 'weight': 4, 'Ao': None, 'Mo': None, 'marked': [], 'lpp': UTCDateTime(2016, 1, 24, 10, 28, 56), 'Mw': None, 'fc': None, 'snr': None, 'epp': UTCDateTime(2016, 1, 24, 10, 28, 24), 'mpp': UTCDateTime(2016, 1, 24, 10, 28, 40), 'w0': None, 'spe': None, 'fm': None, 'channel': u'LHE'}}
        with HidePrints():
            result, station = autopickstation(wfstream=self.a106, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin=self.origin)
        self.assertEqual(expected, result)

if __name__ == '__main__':
    unittest.main()
