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
    def __enter__(self):
        self._original_stdout = sys.stdout
        devnull = open(os.devnull, "w")
        sys.stdout = devnull

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout


class MockParser:

    @staticmethod
    def get_coordinates(station_id):
        """
        Mocks the method get_coordinates from obspy.io.xseed.parser.Parser object
        to avoid building a parser for the unit tests
        :param station_id: 'GR.GRA1..LHZ' or similar
        :type station_id: str
        :return: dictionary containing azimuth, dip, elevation, latitude, longitude,
        local depth as keys
        :rtype: dict

        >>>MockParser.get_coordinates('GR.GRA2')
        {u'azimuth': 0.0, u'dip': -90.0, u'elevation': 512.0, u'latitude': 49.655208, u'local_depth': 0.0, u'longitude': 11.359444}
        """
        station_names = ['GR.GRA1', 'GR.GRA2', 'G.ECH']
        gra1 = {u'azimuth': 0.0, u'dip': -90.0, u'elevation': 499.5, u'latitude': 49.691888, u'local_depth': 0.0, u'longitude': 11.22172}
        gra2 = {u'azimuth': 0.0, u'dip': -90.0, u'elevation': 512.0, u'latitude': 49.655208, u'local_depth': 0.0, u'longitude': 11.359444}
        ech = {u'azimuth': 90.0, u'dip': 0.0, u'elevation': 580.0, u'latitude': 48.216313, u'local_depth': 250.0, u'longitude': 7.158961}
        coordinates = [gra1, gra2, ech]

        for index, name in enumerate(station_names):
            if station_id.startswith(name):
                return coordinates[index]


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
        self.metadata = ('dless', MockParser())

    @staticmethod
    def construct_error_message(station_id, taupy_status, expected, got, phase=''):
        error_message = """Difference in auto {phase} picks for station {station_id} with taupy {status}.\n
                        Expeceted: {exp}\n
                        Got: {got}\n""".format(station_id=station_id, status='enabled' if taupy_status else 'disabled', phase=phase, exp=expected, got=got)
        return error_message

    #@skip("Works")
    def test_autopickstation_taupy_disabled_gra1(self):
        expected = {'P': {'picker': 'auto', 'snrdb': 15.405649120980094, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 34.718816470730317, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 31, 690000), 'w0': None, 'spe': 0.93333333333333235, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 28, 890000), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 32, 690000), 'fm': 'D', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.669661906545489, 'network': u'GR', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 30, 690000), 'snr': 11.667187857573905, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 21, 690000), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 29, 690000), 'fm': None, 'spe': 2.6666666666666665, 'channel': u'LHE'}, 'station': u'GRA1'}
        with HidePrints():
            result = autopickstation(wfstream=self.gra1, pickparam=self.pickparam_taupy_disabled, metadata=(None, None))
        error_message = self.construct_error_message('GR.GRA1', False, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

    #@skip("Works")
    def test_autopickstation_taupy_enabled_gra1(self):
        expected = {'P': {'picker': 'auto', 'snrdb': 15.599905299126778, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 36.307013769185403, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 27, 690000), 'w0': None, 'spe': 0.93333333333333235, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 24, 890000), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 28, 690000), 'fm': 'U', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.669661906545489, 'network': u'GR', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 30, 690000), 'snr': 11.667187857573905, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 21, 690000), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 29, 690000), 'fm': None, 'spe': 2.6666666666666665, 'channel': u'LHE'}, 'station': u'GRA1'}
        with HidePrints():
            result = autopickstation(wfstream=self.gra1, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin=self.origin)
        error_message = self.construct_error_message('GR.GRA1', True, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

    #@skip("works")
    def test_autopickstation_taupy_disabled_gra2(self):
        expected = {'P': {'picker': 'auto', 'snrdb': None, 'weight': 9, 'Mo': None, 'marked': 'shortsignallength', 'Mw': None, 'fc': None, 'snr': None, 'mpp': UTCDateTime(2016, 1, 24, 10, 36, 59, 150000), 'w0': None, 'spe': None, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 36, 43, 150000), 'lpp': UTCDateTime(2016, 1, 24, 10, 37, 15, 150000), 'fm': 'N', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': u'GR', 'weight': 4, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 37, 15, 150000), 'snr': None, 'epp': UTCDateTime(2016, 1, 24, 10, 36, 43, 150000), 'mpp': UTCDateTime(2016, 1, 24, 10, 36, 59, 150000), 'fm': None, 'spe': None, 'channel': u'LHE'}, 'station': u'GRA2'}
        with HidePrints():
            result = autopickstation(wfstream=self.gra2, pickparam=self.pickparam_taupy_disabled, metadata=(None, None))
        error_message = self.construct_error_message('GR.GRA2', False, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

    #skip("")
    def test_autopickstation_taupy_enabled_gra2(self):
        expected = {'P': {'picker': 'auto', 'snrdb': 13.957959025719253, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 24.876879503607871, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 29, 150000), 'w0': None, 'spe': 1.0, 'network': u'GR', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 26, 150000), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 30, 150000), 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 10.573236990555648, 'network': u'GR', 'weight': 1, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 34, 150000), 'snr': 11.410999834108294, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 21, 150000), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 33, 150000), 'fm': None, 'spe': 4.666666666666667, 'channel': u'LHE'}, 'station': u'GRA2'}
        with HidePrints():
            result = autopickstation(wfstream=self.gra2, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin = self.origin)
        error_message = self.construct_error_message('GR.GRA2', True, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

    #@skip("Works")
    def test_autopickstation_taupy_disabled_ech(self):
        expected = {'P': {'picker': 'auto', 'snrdb': None, 'weight': 9, 'Mo': None, 'marked': 'SinsteadP', 'Mw': None, 'fc': None, 'snr': None, 'mpp': UTCDateTime(2016, 1, 24, 10, 26, 57), 'w0': None, 'spe': None, 'network': u'G', 'epp': UTCDateTime(2016, 1, 24, 10, 26, 41), 'lpp': UTCDateTime(2016, 1, 24, 10, 27, 13), 'fm': 'N', 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': None, 'network': u'G', 'weight': 4, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 27, 13), 'snr': None, 'epp': UTCDateTime(2016, 1, 24, 10, 26, 41), 'mpp': UTCDateTime(2016, 1, 24, 10, 26, 57), 'fm': None, 'spe': None, 'channel': u'LHE'}, 'station': u'ECH'}
        with HidePrints():
            result = autopickstation(wfstream=self.ech, pickparam=self.pickparam_taupy_disabled)
        error_message = self.construct_error_message('G.ECH', False, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

    #@skip("Works")
    def test_autopickstation_taupy_enabled_ech(self):
        # this station has a long time of before the first onset, so taupy will help during picking
        expected = {'P': {'picker': 'auto', 'snrdb': 9.9753586609166316, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 9.9434218804137107, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 34), 'w0': None, 'spe': 1.6666666666666667, 'network': u'G', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 29), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 35), 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 12.698999454169567, 'network': u'G', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 44), 'snr': 18.616581906366577, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 33), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 43), 'fm': None, 'spe': 3.3333333333333335, 'channel': u'LHE'}, 'station': u'ECH'}
        with HidePrints():
            result = autopickstation(wfstream=self.ech, pickparam=self.pickparam_taupy_enabled, metadata=self.metadata, origin=self.origin)
        error_message = self.construct_error_message('G.ECH', True, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

    @skip("Results are just a copy of ech")
    def test_autopickstation_taupy_disabled_fiesa(self):
        # this station has a long time of before the first onset, so taupy will help during picking
        expected = {'P': {'picker': 'auto', 'snrdb': 9.9753586609166316, 'weight': 0, 'Mo': None, 'marked': [], 'Mw': None, 'fc': None, 'snr': 9.9434218804137107, 'mpp': UTCDateTime(2016, 1, 24, 10, 41, 34), 'w0': None, 'spe': 1.6666666666666667, 'network': u'G', 'epp': UTCDateTime(2016, 1, 24, 10, 41, 29), 'lpp': UTCDateTime(2016, 1, 24, 10, 41, 35), 'fm': None, 'channel': u'LHZ'}, 'S': {'picker': 'auto', 'snrdb': 12.698999454169567, 'network': u'G', 'weight': 0, 'Ao': None, 'lpp': UTCDateTime(2016, 1, 24, 10, 50, 44), 'snr': 18.616581906366577, 'epp': UTCDateTime(2016, 1, 24, 10, 50, 33), 'mpp': UTCDateTime(2016, 1, 24, 10, 50, 43), 'fm': None, 'spe': 3.3333333333333335, 'channel': u'LHE'}, 'station': u'ECH'}
        with HidePrints():
            result = autopickstation(wfstream=self.fiesa, pickparam=self.pickparam_taupy_disabled)
        error_message = self.construct_error_message('CH.FIESE', True, expected, result)
        self.assertDictContainsSubset(expected=expected['P'], actual=result['P'])
        self.assertDictContainsSubset(expected=expected['S'], actual=result['S'])
        self.assertEqual(expected['station'], result['station'])

if __name__ == '__main__':
    unittest.main()
