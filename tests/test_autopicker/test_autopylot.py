import os
import pytest

from obspy import read_events

from autoPyLoT import autoPyLoT


class TestAutopickerGlobal():
    def init(self):
        self.params_infile = 'pylot_alparray_mantle_corr_stack_0.03-0.5.in'
        self.test_event_dir = 'dmt_database_test'
        self.fname_outfile_xml = os.path.join(
            self.test_event_dir, '20171010_063224.a', 'PyLoT_20171010_063224.a_autopylot.xml'
        )

        # check if the input files exist
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
        # check for output file in test directory and remove it if necessary
        if os.path.isfile(self.fname_outfile_xml):
            os.remove(self.fname_outfile_xml)
        autoPyLoT(inputfile=self.params_infile, eventid='20171010_063224.a', obspyDMT_wfpath='processed')

        # test for different known output files if they are identical or not
        compare_pickfiles(self.fname_outfile_xml, 'PyLoT_20171010_063224.a_autopylot.xml', True)
        compare_pickfiles(self.fname_outfile_xml, 'PyLoT_20171010_063224.a_saved_from_GUI.xml', True)
        compare_pickfiles(self.fname_outfile_xml, 'PyLoT_20171010_063224.a_corrected_taup_times_0.03-0.5_P.xml', False)


def compare_pickfiles(pickfile1: str, pickfile2: str, samefile: bool = True) -> None:
    """
    Compare the pick times and errors from two pick files.

    Parameters:
        pickfile1 (str): The path to the first pick file.
        pickfile2 (str): The path to the second pick file.
        samefile (bool): A flag indicating whether the two files are expected to be the same. Defaults to True.

    Returns:
        None
    """
    cat1 = read_events(pickfile1)
    cat2 = read_events(pickfile2)
    picks1 = sorted(cat1[0].picks, key=lambda pick: str(pick.waveform_id))
    picks2 = sorted(cat2[0].picks, key=lambda pick: str(pick.waveform_id))
    pick_times1 = [pick.time for pick in picks1]
    pick_times2 = [pick.time for pick in picks2]
    pick_terrs1 = [pick.time_errors for pick in picks1]
    pick_terrs2 = [pick.time_errors for pick in picks2]

    # check if times and errors are identical or not depending on the samefile flag
    assert (pick_times1 == pick_times2) is samefile, 'Pick times error'
    assert (pick_terrs1 == pick_terrs2) is samefile, 'Pick time errors errors'
