import pytest
from obspy import read, Trace, UTCDateTime

from pylot.correlation.pick_correlation_correction import XCorrPickCorrection


class TestXCorrPickCorrection():
    def setup(self):
        self.make_test_traces()
        self.make_test_picks()
        self.t_before = 2.
        self.t_after = 2.
        self.cc_maxlag = 0.5

    def make_test_traces(self):
        # take first trace of test Stream from obspy
        tr1 = read()[0]
        # filter trace
        tr1.filter('bandpass', freqmin=1, freqmax=20)
        # make a copy and shift the copy by 0.1 s
        tr2 = tr1.copy()
        tr2.stats.starttime += 0.1

        self.trace1 = tr1
        self.trace2 = tr2

    def make_test_picks(self):
        # create an artificial reference pick on reference trace (trace1) and another one on the 0.1 s shifted trace
        self.tpick1 = UTCDateTime('2009-08-24T00:20:07.7')
        # shift the second pick by 0.2 s, the correction should be around 0.1 s now
        self.tpick2 = self.tpick1 + 0.2

    def test_slice_trace_okay(self):

        self.setup()
        xcpc = XCorrPickCorrection(UTCDateTime(), Trace(), UTCDateTime(), Trace(),
                                   t_before=self.t_before, t_after=self.t_after, cc_maxlag=self.cc_maxlag)

        test_trace = self.trace1
        pick_time = self.tpick2

        sliced_trace = xcpc.slice_trace(test_trace, pick_time)
        assert ((sliced_trace.stats.starttime == pick_time - self.t_before - self.cc_maxlag / 2)
                and (sliced_trace.stats.endtime == pick_time + self.t_after + self.cc_maxlag / 2))

    def test_slice_trace_fails(self):
        self.setup()

        test_trace = self.trace1
        pick_time = self.tpick1

        with pytest.raises(Exception):
            xcpc = XCorrPickCorrection(UTCDateTime(), Trace(), UTCDateTime(), Trace(),
                                       t_before=self.t_before - 20, t_after=self.t_after, cc_maxlag=self.cc_maxlag)
            xcpc.slice_trace(test_trace, pick_time)

        with pytest.raises(Exception):
            xcpc = XCorrPickCorrection(UTCDateTime(), Trace(), UTCDateTime(), Trace(),
                                       t_before=self.t_before, t_after=self.t_after + 50, cc_maxlag=self.cc_maxlag)
            xcpc.slice_trace(test_trace, pick_time)

    def test_cross_correlation(self):
        self.setup()

        # create XCorrPickCorrection object
        xcpc = XCorrPickCorrection(self.tpick1, self.trace1, self.tpick2, self.trace2, t_before=self.t_before,
                                   t_after=self.t_after, cc_maxlag=self.cc_maxlag)

        # execute correlation
        correction, cc_max, uncert, fwfm = xcpc.cross_correlation(False, '', '')

        # define awaited test result
        test_result = (-0.09983091718314982, 0.9578431835689154, 0.0015285160561610929, 0.03625786256084631)

        # check results
        assert pytest.approx(test_result, rel=1e-6) == (correction, cc_max, uncert, fwfm)