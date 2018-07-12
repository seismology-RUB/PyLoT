import unittest
import os

from obspy import UTCDateTime
from pylot.core.util.dataprocessing import Metadata

class TestMetadata(unittest.TestCase):

    def setUp(self):
        self.station_id = 'BW.WETR..HH'
        self.time = UTCDateTime('2012-08-01')
        metadata_folder = 'metadata1'
        self.m = Metadata(metadata_folder)

    def test_get_coordinates_sucess(self):
        expected = {'Z': {u'elevation': 607.0, u'longitude': 12.87571, u'local_depth': 0.0, u'azimuth': 0.0, u'latitude': 49.14502, u'dip': -90.0},
                    'E': {u'azimuth': 90.0, u'dip': 0.0, u'elevation': 607.0, u'latitude': 49.14502, u'local_depth': 0.0, u'longitude': 12.87571},
                    'N': {u'azimuth': 0.0, u'dip': 0.0, u'elevation': 607.0, u'latitude': 49.14502, u'local_depth': 0.0, u'longitude': 12.87571}
                    }
        result = {}
        for channel in ('Z', 'N', 'E'):
            coords = self.m.get_coordinates(self.station_id+channel, time=self.time)
            result[channel] = coords
            self.assertDictEqual(result[channel], expected[channel])

class TestMetadataAdding(unittest.TestCase):
    """Tests if adding files and directories to a metadata object works."""

    def setUp(self):
        self.station_id = 'BW.WETR..HH'
        self.metadata_folders = ('metadata1', 'metadata2')
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
        self.assertEqual(['metadata1/DATALESS.BW.WETR..HHZ'], self.m.inventory_files.keys())  # does the filename exist in inventory files?
        self.assertEqual(['data', 'invtype'], self.m.inventory_files['metadata1/DATALESS.BW.WETR..HHZ'].keys())  # is the required information attacht to the filename?
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
        self.metadata_folders = ('metadata1', 'metadata2')
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
        self.m.remove_inventory('metadata_not_existing')
        self.assertIn(self.metadata_folders[0], self.m.inventories)

    def isEmpty(self, metadata):
        """Asserts if the given metadata object is empty"""
        self.assertDictEqual({}, metadata.inventory_files)
        self.assertDictEqual({}, metadata.seed_ids)
        self.assertEqual([], metadata.inventories)

