#!/usr/bin/env python3

import pandas as pd
import time
import os
import argparse
import traceback
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoServiceException

"""
seismo_fetch.py

Query seismic data and station metadata saved to files.

Usage (CLI):
    python seismo_fetch.py --fname quries.csv --client IRIS --outdir ./data --format mseed

Import:
    from seismo_fetch import all

Last modified by spipatprathanporn@ucsd.edu, 09/19/2026
"""

def _cli():
    p = argparse.ArgumentParser(description="Query seismic data and station metadata saved to files.")
    p.add_argument("--fname", type=str, required=True, help="File name of the query table")
    p.add_argument("--client", type=str, default="IRIS", help="FDSN client name (IRIS, GFZ, etc.)")
    p.add_argument("--outdir", type=str, default=".", help="output directory")
    p.add_argument("--format", type=str, default="mseed", choices=["mseed", "sac"])
    args = p.parse_args()
    
    client = Client(args.client)
    outdir = args.outdir
    format = args.format

    df = pd.read_csv(args.fname, dtype={"location": str})

    os.makedirs(outdir, exist_ok=True)

    inv = {}
    for index, row in df.iterrows():
        netcode = row['network']
        stacode = row['station']
        loc = row['location']
        chan = row['channel']
        starttime = UTCDateTime(row['starttime'])
        endtime = UTCDateTime(row['endtime'])

        # station geographical position
        staname = f"{netcode}.{stacode}"
        # find station if it does not exist in the inventory yet
        if not staname in inv:
            dat = client.get_stations(starttime=starttime, endtime=endtime,
                                      network=netcode, station=stacode, level="station")
            inv[staname] = {"latitude": dat[0][0].latitude, 
                            "longitude": dat[0][0].longitude, 
                            "elevation": dat[0][0].elevation}
        
        retries = 3
        delay = 1
        skip_index = False
        for attempt in range(retries):
            print(f"Downloading {netcode}.{stacode}.{loc}.{chan}...")
            # try to download and save SACPZ response for this channel (if available)
            try:
                inv_resp = client.get_stations(network=netcode, station=stacode,
                                                location=loc or "*", channel=chan,
                                                starttime=starttime, endtime=endtime,
                                                level="response")
                resp_fname = f"{netcode}.{stacode}.{loc}.{chan}.{starttime.isoformat()}_{endtime.isoformat()}.sacpz"
                resp_path = os.path.join(outdir, resp_fname)
                inv_resp.write(resp_path, format="SACPZ")
                print("response downloaded:", resp_path)
                print("RESPONSE SAVED:", resp_path)
                # Be a polite scraper: rest 0.5 seconds between iterations
                time.sleep(0.5)
                break
            except FDSNNoServiceException or ValueError:
                if attempt < retries - 1:
                    print(f"Throttled. Retrying row {index} in {delay}s...")
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                    continue
                print("failed to download response for:", netcode, stacode, loc, chan)
            except Exception:
                traceback.print_exc()
                print("no response available for:", netcode, stacode, loc, chan)
                skip_index = True
                break
            # a successful run should exit early in the TRY block
            skip_index = True
            break
        if skip_index:
             continue

        retries = 3
        delay = 1
        skip_index = False
        for attempt in range(retries):
            try:
                tr = client.get_waveforms(netcode, stacode, loc, chan, starttime, endtime)
                
                # Be a polite scraper: rest 0.5 seconds between iterations
                time.sleep(0.5)
                break
            except FDSNNoServiceException or ValueError:
                if attempt < retries - 1:
                    print(f"Throttled. Retrying row {index} in {delay}s...")
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                    continue
                print("failed to download waveform for:", netcode, stacode, loc, chan)
            except Exception:
                    traceback.print_exc()
                    print("no waveform available for:", netcode, stacode, loc, chan)
                    skip_index = True
                    break
            skip_index = True
            break
        if skip_index:
             continue

        if len(tr) <= 0:
            print("failed to download waveform for:", netcode, stacode, loc, chan)
            continue
        else:
            print("waveform downloaded:", netcode, stacode, loc, chan)

        # write SAC (first trace) and populate SAC header with station metadata
        if format.lower() == "sac":
            fname = f"{netcode}.{stacode}.{loc}.{chan}.{starttime.isoformat()}_{endtime.isoformat()}.sac"
            path = os.path.join(outdir, fname)
            tr0 = tr[0]
            # ensure SAC header dict exists
            if not hasattr(tr0.stats, "sac") or tr0.stats.sac is None:
                    tr0.stats.sac = {}
            sac = tr0.stats.sac
            # station metadata
            sac["stla"] = inv[staname]["latitude"]
            sac["stlo"] = inv[staname]["longitude"]
            sac["stel"] = inv[staname]["elevation"]
            sac["knetwk"] = netcode
            sac["kstnm"] = stacode
            sac["kcmpnm"] = chan
            sac["khole"] = loc or ""
            tr0.write(path, format="SAC")
            print("SAVED:", path)
        else:
            # mseed: Combine stream into one file per channel/time
            fname = f"{netcode}.{stacode}.{loc}.{chan}.{starttime.isoformat()}_{endtime.isoformat()}.mseed"
            path = os.path.join(outdir, fname)
            tr.write(path, format="MSEED")
            print("SAVED:", path)

if __name__ == "__main__":
    _cli()