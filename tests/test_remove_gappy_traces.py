import os
import obspy
DATADIR = os.path.join(os.path.dirname(__file__), 'waveforms')

def test_can_handle_data_with_duplicate_packets():
    from wphase.data_acquisition import remove_gappy_traces
    st0 = obspy.read(DATADIR + '/duplicates.mseed')
    st1 = remove_gappy_traces(st0)
    assert len(set(tr.meta.station for tr in st1)) == 4
