"""CLI entrypoint to fetch inventory for an event from FDSN.

Usage: python3 -m wphase.fetch_inventory <server> <evid> <output_path>"""

import sys
from obspy.clients.fdsn import Client
from wphase._runner_fdsn import load_metadata
from wphase.psi.model import Event

def fetch_inventory(
    server: str,
    evid: str,
    save_path: str,
    dist_range=[5,90],
    networks="ALL",
):
    client = Client(server)
    cat = client.get_events(eventid=evid)
    origin = cat.events[0].preferred_origin()
    eqinfo = Event(
        longitude=origin.longitude,
        latitude=origin.latitude,
        depth=origin.depth / 1000,
        time=origin.time,
    )
    load_metadata(
        client=client,
        eqinfo=eqinfo,
        dist_range=dist_range,
        networks=networks,
        save_path=save_path, 
    )

if __name__ == "__main__":
    fetch_inventory(*sys.argv[1:])
