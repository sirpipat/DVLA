#!/usr/bin/env python3
"""
seismo_query.py

Query seismic stations within `radius` degrees of (lat, lon) that were active during
[starttime, endtime] and optionally download waveforms.

Usage (CLI):
    python seismo_query.py --lat 34.05 --lon -118.25 --radius 2 \
            --start 2020-01-01T00:00:00 --end 2020-01-01T01:00:00 \
            --client IRIS --channels BH? --outdir ./data --format mseed

Import:
    from seismo_query import find_stations, download_waveforms, parse_time
"""

from typing import List, Tuple
import os
import argparse
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.inventory import Inventory

def parse_time(t: str) -> UTCDateTime:
        return UTCDateTime(t)

def find_stations(client: Client, lat: float, lon: float, radius_deg: float,
                                    starttime: UTCDateTime, endtime: UTCDateTime,
                                    network: str = "*", location: str = "*", channel: str = "*",
                                    level: str = "channel") -> Inventory:
        """
        Return an Inventory (StationXML) of stations within `radius_deg` of (lat,lon)
        that have availability intersecting [starttime, endtime].
        """
        return client.get_stations(latitude=lat, longitude=lon, maxradius=radius_deg,
                                                             starttime=starttime, endtime=endtime,
                                                             network=network, location=location, channel=channel,
                                                             level=level)

def download_waveforms(client: Client, inventory: Inventory,
                                             starttime: UTCDateTime, endtime: UTCDateTime,
                                             outdir: str = ".", fmt: str = "mseed") -> List[str]:
        """
        Download waveforms for all channels in `inventory` between starttime and endtime.
        Saves files to outdir. Returns list of saved file paths.
        fmt: "mseed" or "sac"
        """

        os.makedirs(outdir, exist_ok=True)
        saved = []
        saved_responses = []
        for net in inventory:
                netcode = net.code
                for sta in net:
                        stacode = sta.code
                        for cha in sta.channels:
                                chan = cha.code
                                loc = cha.location_code or ""
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
                                        saved_responses.append(resp_path)
                                        print("response downloaded:", resp_path)
                                except Exception:
                                        print("no response available for:", netcode, stacode, loc, chan)

                                try:
                                        tr = client.get_waveforms(netcode, stacode, loc, chan, starttime, endtime)
                                        print("waveform downloaded:", netcode, stacode, loc, chan)
                                except Exception:
                                        print("failed to download:", netcode, stacode, loc, chan)
                                        continue
                                # write SAC (first trace) and populate SAC header with station metadata
                                if fmt.lower() == "sac":
                                        fname = f"{netcode}.{stacode}.{loc}.{chan}.{starttime.isoformat()}_{endtime.isoformat()}.sac"
                                        path = os.path.join(outdir, fname)
                                        tr0 = tr[0]
                                        # ensure SAC header dict exists
                                        if not hasattr(tr0.stats, "sac") or tr0.stats.sac is None:
                                                tr0.stats.sac = {}
                                        sac = tr0.stats.sac
                                        # station metadata
                                        sac["stla"] = getattr(sta, "latitude", 0.0)
                                        sac["stlo"] = getattr(sta, "longitude", 0.0)
                                        sac["stel"] = getattr(sta, "elevation", 0.0)
                                        sac["knetwk"] = netcode
                                        sac["kstnm"] = sta.code
                                        sac["kcmpnm"] = chan
                                        sac["khole"] = loc or ""
                                        tr0.write(path, format="SAC")
                                        saved.append(path)
                                else:
                                        # mseed: Combine stream into one file per channel/time
                                        fname = f"{netcode}.{stacode}.{loc}.{chan}.{starttime.isoformat()}_{endtime.isoformat()}.mseed"
                                        path = os.path.join(outdir, fname)
                                        tr.write(path, format="MSEED")
                                        saved.append(path)
        # report saved response files
        for r in saved_responses:
                print("RESPONSE SAVED:", r)
        return saved

def _cli():
        p = argparse.ArgumentParser(description="Query stations within radius and download waveforms.")
        p.add_argument("--lat", type=float, required=True)
        p.add_argument("--lon", type=float, required=True)
        p.add_argument("--radius", type=float, required=True, help="degrees")
        p.add_argument("--start", type=str, required=True, help="ISO time")
        p.add_argument("--end", type=str, required=True, help="ISO time")
        p.add_argument("--client", type=str, default="IRIS", help="FDSN client name (IRIS, GFZ, etc.)")
        p.add_argument("--network", type=str, default="*", help="network code or *")
        p.add_argument("--location", type=str, default="*", help="location code or *")
        p.add_argument("--channels", type=str, default="BH?", help="channel code or pattern (e.g. BH?)")
        p.add_argument("--outdir", type=str, default=".", help="output directory")
        p.add_argument("--format", type=str, default="mseed", choices=["mseed", "sac"])
        args = p.parse_args()

        client = Client(args.client)
        start = parse_time(args.start)
        end = parse_time(args.end)
        inv = find_stations(client, args.lat, args.lon, args.radius, start, end,
                                                network=args.network, location=args.location, channel=args.channels)
        # print a concise summary
        for net in inv:
                for sta in net:
                        print(f"{net.code}.{sta.code}  ({sta.latitude},{sta.longitude})")
        saved = download_waveforms(client, inv, start, end, outdir=args.outdir, fmt=args.format)
        for s in saved:
                print("SAVED:", s)

if __name__ == "__main__":
        _cli()