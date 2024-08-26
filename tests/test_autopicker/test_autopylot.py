import os
import pytest

from autoPyLoT import autoPyLoT


class TestAutopickerGlobal():
    def init(self):
        self.params_infile = 'pylot_alparray_mantle_corr_stack_0.03-0.5.in'
        self.test_event_dir = 'dmt_database_test'

        if not os.path.isfile(self.params_infile):
            print(f'Test input file {os.path.abspath(self.params_infile)} not found.')
            return False

        if not os.path.exists(self.test_event_dir):
            print(
                f'Test event directory not found at location "{os.path.abspath(self.test_event_dir)}". '
                f'Make sure to load it from the website first.'
            )
            return False

        return True

    def test_autopicker(self):
        assert self.init(), 'Initialization failed due to missing input files.'
        #autoPyLoT(inputfile=self.params_infile, eventid='20171010_063224.a')
