import unittest
import os

from obspy import UTCDateTime
from obspy.io.xseed.utils import SEEDParserException
from obspy.io.xseed import Parser
from pylot.core.util.dataprocessing import Metadata
from tests.utils import HidePrints

class TestMetadata(unittest.TestCase):

    def setUp(self):
        self.station_id = 'BW.WETR..HH'
        self.time = UTCDateTime('2012-08-01')
        metadata_folder = os.path.join('test_data', 'dless_multiple_files', 'metadata1')
        self.m = Metadata(metadata_folder)

    def test_get_coordinates_sucess(self):
        expected = {'Z': {u'elevation': 607.0, u'longitude': 12.87571, u'local_depth': 0.0, u'azimuth': 0.0, u'latitude': 49.14502, u'dip': -90.0},
                    'E': {u'azimuth': 90.0, u'dip': 0.0, u'elevation': 607.0, u'latitude': 49.14502, u'local_depth': 0.0, u'longitude': 12.87571},
                    'N': {u'azimuth': 0.0, u'dip': 0.0, u'elevation': 607.0, u'latitude': 49.14502, u'local_depth': 0.0, u'longitude': 12.87571}
                    }
        result = {}
        for channel in ('Z', 'N', 'E'):
            with HidePrints():
                coords = self.m.get_coordinates(self.station_id+channel, time=self.time)
            result[channel] = coords
            self.assertDictEqual(result[channel], expected[channel])

    def test_get_coordinates_sucess_no_time(self):
        expected = {'Z': {u'elevation': 607.0, u'longitude': 12.87571, u'local_depth': 0.0, u'azimuth': 0.0, u'latitude': 49.14502, u'dip': -90.0},
                    'E': {u'azimuth': 90.0, u'dip': 0.0, u'elevation': 607.0, u'latitude': 49.14502, u'local_depth': 0.0, u'longitude': 12.87571},
                    'N': {u'azimuth': 0.0, u'dip': 0.0, u'elevation': 607.0, u'latitude': 49.14502, u'local_depth': 0.0, u'longitude': 12.87571}
                    }
        result = {}
        for channel in ('Z', 'N', 'E'):
            with HidePrints():
                coords = self.m.get_coordinates(self.station_id+channel)
            result[channel] = coords
            self.assertDictEqual(result[channel], expected[channel])


class TestMetadataAdding(unittest.TestCase):
    """Tests if adding files and directories to a metadata object works."""

    def setUp(self):
        self.station_id = 'BW.WETR..HH'
        self.metadata_folders = (os.path.join('test_data', 'dless_multiple_files', 'metadata1'),
                                 os.path.join('test_data', 'dless_multiple_files', 'metadata2'))
        self.m = Metadata()

    def test_add_inventory_folder(self):
        """Test if add_inventory adds the folder to the list of inventories"""
        self.m.add_inventory(self.metadata_folders[0])
        # adding an inventory folder should append it to the list of inventories
        self.assertDictEqual({}, self.m.inventory_files)
        self.assertDictEqual({}, self.m.seed_ids)
        self.assertEqual([self.metadata_folders[0]], self.m.inventories)

    def test_add_inventory_file(self):
        """Test if add_inventory_file adds the folder containing the file to the list of inventories and
        if the files is added to inventory_files"""
        fpath = os.path.join(self.metadata_folders[0], 'DATALESS.BW.WETR..HHZ')
        self.m.add_inventory_file(fpath)
        # adding an inventory file should append its folder to the list of inventories and the file to the
        self.assertEqual([os.path.join(self.metadata_folders[0], 'DATALESS.BW.WETR..HHZ')], self.m.inventory_files.keys())  # does the filename exist in inventory files?
        self.assertEqual(['data', 'invtype'], self.m.inventory_files[os.path.join(self.metadata_folders[0], 'DATALESS.BW.WETR..HHZ')].keys())  # is the required information attacht to the filename?
        self.assertDictEqual({}, self.m.seed_ids)
        self.assertEqual([self.metadata_folders[0]], self.m.inventories)

    def test_add_inventory_invalid_path(self):
        """Test if adding an inventory that is not an existing directory fails with an exception"""
        with self.assertRaises(Exception):
            self.m.add_inventory('InvalidDirName')
        self.assertEqual([], self.m.inventories)  # inventory list should still be empty

    def test_add_inventory_file_invalid_path(self):
        """Test if adding a inventory file with an invalid path fails with an exception"""
        with self.assertRaises(Exception):
            self.m.add_inventory_file('/invalid/file/name')
        self.assertEqual([], self.m.inventories)  # inventory list should still be empty


class TestMetadataRemoval(unittest.TestCase):
    """Tests if removing files and directories to a metadata object works."""

    def setUp(self):
        self.station_id = 'BW.WETR..HH'
        self.metadata_folders = (os.path.join('test_data', 'dless_multiple_files', 'metadata1'),
                                 os.path.join('test_data', 'dless_multiple_files', 'metadata2'))
        self.m = Metadata()

    def test_remove_all_inventories(self):
        """Test if function remove_inventory cleans the Metadata object """
        # add multiple inventories
        for folder in self.metadata_folders:
            self.m.add_inventory(folder)
        self.m.remove_all_inventories()
        self.isEmpty(self.m)

    def test_remove_inventory(self):
        """Test if remove_inventory removes single inventories"""
        # add multiple inventories
        for folder in self.metadata_folders:
            self.m.add_inventory(folder)
        self.m.remove_inventory(self.metadata_folders[0])
        self.assertNotIn(self.metadata_folders[0], self.m.inventories)
        self.m.remove_inventory(self.metadata_folders[1])
        self.assertNotIn(self.metadata_folders[1], self.m.inventories)
        self.isEmpty(self.m)

    def test_remove_inventory_not_in_inventory_list(self):
        """Test if remove_inventory does not modify the metadata instance if the given inventory to remove does not
        exist in the instance."""
        # add multiple inventories
        self.m.add_inventory(self.metadata_folders[0])
        with HidePrints():
            self.m.remove_inventory('metadata_not_existing')
        self.assertIn(self.metadata_folders[0], self.m.inventories)

    def isEmpty(self, metadata):
        """Asserts if the given metadata object is empty"""
        self.assertDictEqual({}, metadata.inventory_files)
        self.assertDictEqual({}, metadata.seed_ids)
        self.assertEqual([], metadata.inventories)


class TestMetadata_read_single_file(unittest.TestCase):

    def setUp(self):
        self.station_id = 'BW.WETR..HHZ'
        self.metadata_folders = (os.path.join('test_data', 'dless_multiple_files', 'metadata1'),
                                 os.path.join('test_data', 'dless_multiple_files', 'metadata2'))
        self.metadata_paths = []
        self.m = Metadata()

    def test_read_single_file(self):
        """Test if reading a single file works"""
        fname = os.path.join(self.metadata_folders[0], 'DATALESS.'+self.station_id)
        with HidePrints():
            res = self.m.read_single_file(fname)
        # method should return true if file is successfully read
        self.assertTrue(res)
        # list of inventories (folders) should be empty
        self.assertEqual([], self.m.inventories)
        # list of inventory files should contain the added file
        self.assertIn(fname, self.m.inventory_files.keys())
        self.assertEqual({}, self.m.seed_ids)

    def test_read_single_file_invalid_path(self):
        """Test if reading from a non existing file fails. The filename should not be
        added to the metadata object"""
        fname = os.path.join("this", "path", "doesnt", "exist")
        with HidePrints():
            res = self.m.read_single_file(fname)
        # method should return None if file reading fails
        self.assertIsNone(res)
        # list of inventories (folders) should be empty
        self.assertEqual([], self.m.inventories)
        # list of inventory files should not contain the added file
        self.assertNotIn(fname, self.m.inventory_files.keys())
        self.assertEqual({}, self.m.seed_ids)

    def test_read_single_file_multiple_times(self):
        """Test if reading a file twice doesnt add it twice to the metadata object"""
        fname = os.path.join(self.metadata_folders[0], 'DATALESS.'+self.station_id)
        with HidePrints():
            res1 = self.m.read_single_file(fname)
            res2 = self.m.read_single_file(fname)
        self.assertTrue(res1)
        self.assertIsNone(res2)
        self.assertItemsEqual([fname], self.m.inventory_files.keys())


class TestMetadataMultipleTime(unittest.TestCase):
    """Test if stations with multiple metadata entries in a single file are handled correctly.
    The user must specify the time where he wants to get metadata.

    The station ROTT changed has metadata available at multiple times
    LE.ROTT..HNE | 200.00 Hz | Titan 4g-EDR-209, Very Low gain, 200 sps | 2015-01-08 - 2015-03-19 | Lat: 49.1, Lng: 8.1
	LE.ROTT..HNE | 200.00 Hz | Titan 4g-EDR-209, Very Low gain, 200 sps | 2015-03-19 -  | Lat: 49.1, Lng: 8.1
	LE.ROTT..HNN | 200.00 Hz | Titan 4g-EDR-209, Very Low gain, 200 sps | 2015-01-08 - 2015-03-19 | Lat: 49.1, Lng: 8.1
	LE.ROTT..HNN | 200.00 Hz | Titan 4g-EDR-209, Very Low gain, 200 sps | 2015-03-19 -  | Lat: 49.1, Lng: 8.1
	LE.ROTT..HNZ | 200.00 Hz | Titan 4g-EDR-209, Very Low gain, 200 sps | 2015-01-08 - 2015-03-19 | Lat: 49.1, Lng: 8.1
	LE.ROTT..HNZ | 200.00 Hz | Titan 4g-EDR-209, Very Low gain, 200 sps | 2015-03-19 -  | Lat: 49.1, Lng: 8.1
    """

    def setUp(self):
        self.seed_id = 'LE.ROTT..HN'
        path = os.path.dirname(__file__)  # gets path to currently running script
        metadata = os.path.join('test_data', 'dless_multiple_times', 'MAGS2_LE_ROTT.dless')  # specific subfolder of test data
        metadata_path = os.path.join(path, metadata)
        self.m = Metadata(metadata_path)
        self.p = Parser(metadata_path)

    def test_get_metadata_works_without_datetime(self):
        """Test if get_metadata works if multiple metadata entries are available but no time is
        specified."""
        for channel in ('Z', 'N', 'E'):
            with HidePrints():
                md = self.m.get_metadata(self.seed_id + channel)
            self.assertDictEqual(md['data'].get_inventory(), self.p.get_inventory())

    def test_get_metadata_works_with_first_datetime(self):
        """Test if get_metadata works if multiple metadata entries are available and the older time is specified."""
        t = UTCDateTime('2015-02-08')
        for channel in ('Z', 'N', 'E'):
            with HidePrints():
                md = self.m.get_metadata(self.seed_id + channel, t)
            self.assertDictEqual(md['data'].get_inventory(), self.p.get_inventory())

    def test_get_metadata_fails_when_time_before_starttime(self):
        """Tests if get_metadata returns None when given a data that is before the start date
        of the metadata"""
        with HidePrints():
            md = self.m.get_metadata(self.seed_id, UTCDateTime('1960-07-20'))
        self.assertIs(md, None)

    def test_get_metadata_invalid_seed_id(self):
        """Tes if get metadata returns none when asked for a seed id that does not exist"""
        with HidePrints():
            res = self.m.get_metadata("this.doesnt..exist")
        self.assertIsNone(res)

        self.assertIsNone(res)
