from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import numpy as np
from scipy.io import savemat
import argparse
import sys
import signal

#!/usr/bin/python3
"""
    Query earthquake catalog with ObsPy and save results to a MATLAB .mat file. It automatically quits if the query takes longer than 10 seconds.
    If you want to query longer time periods, consider splitting your queries into smaller time windows or query manually in python.

    Usage:
        python irisfetch_events.py --starttime '2020-01-01T00:00:00' --endtime 2020-01-02T00:00:00 --minmagnitude 5.0 --output eq_data.mat

    Last modified by spipatprathanporn@ucsd.edu, 2025-12-09
"""

def timeout_handler(signum, frame):
    raise KeyboardInterrupt("Query timed out")

def query_and_save(
        starttime=None,
        endtime=None,
        minlatitude=None,
        maxlatitude=None,
        minlongitude=None,
        maxlongitude=None,
        latitude=None,
        longitude=None,
        minradius=None,
        maxradius=None,
        mindepth=None,
        maxdepth=None,
        minmagnitude=None,
        maxmagnitude=None,
        eventid=None,
        output='events.mat'
):
    try:
        # create FDSN client
        client = Client("IRIS")

        # query events
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(10)  # set timeout to 10 seconds

        catalog = client.get_events(
            starttime=UTCDateTime(starttime) if starttime else None,
            endtime=UTCDateTime(endtime) if endtime else None,
            minlatitude=minlatitude,
            maxlatitude=maxlatitude,
            minlongitude=minlongitude,
            maxlongitude=maxlongitude,
            latitude=latitude,
            longitude=longitude,
            minradius=minradius,
            maxradius=maxradius,
            mindepth=mindepth,
            maxdepth=maxdepth,
            minmagnitude=minmagnitude,
            maxmagnitude=maxmagnitude,
            eventid=eventid
        )
        signal.alarm(0)  # disable alarm
    except  KeyboardInterrupt:
        print("Query timed out after 10 seconds", file=sys.stderr)
        sys.exit(2)
    except Exception as e:
        signal.alarm(0)  # disable alarm in case of exception
        print(f"Error querying events: {e}", file=sys.stderr)
        sys.exit(1)
    
    # process events
    events_data = []
    for event in catalog:
        event_info = process_event(event)
        events_data.append(event_info)

    # save to .mat file
    savemat(output, {'events': events_data})

"""
    Parse command-line arguments.
"""
def parse_args():
    parser = argparse.ArgumentParser(description="Query earthquake catalog with ObsPy and save results to a MATLAB .mat file.")
    parser.add_argument('--starttime', help='Start time in ISO format, e.g., 2020-01-01T00:00:00')
    parser.add_argument('--endtime', help='End time in ISO format, e.g., 2020-01-02T00:00:00')
    parser.add_argument('--minlatitude', type=float, help='Minimum latitude')
    parser.add_argument('--maxlatitude', type=float, help='Maximum latitude')
    parser.add_argument('--minlongitude', type=float, help='Minimum longitude')
    parser.add_argument('--maxlongitude', type=float, help='Maximum longitude')
    parser.add_argument('--latitude', type=float, help='Latitude for circular search')
    parser.add_argument('--longitude', type=float, help='Longitude for circular search')
    parser.add_argument('--minradius', type=float, help='Minimum radius for circular search')
    parser.add_argument('--maxradius', type=float, help='Maximum radius for circular search')
    parser.add_argument('--mindepth', type=float, help='Minimum depth in km')
    parser.add_argument('--maxdepth', type=float, help='Maximum depth in km')
    parser.add_argument('--minmagnitude', type=float, help='Minimum magnitude')
    parser.add_argument('--maxmagnitude', type=float, help='Maximum magnitude')
    parser.add_argument('--eventid', help='Event ID')
    parser.add_argument('--output', required=True, help='Output .mat file name')
    return parser.parse_args()

"""
    Process a single event and extract relevant information.
"""
def process_event(event):
    origin = event.preferred_origin() or event.origins[0]
    magnitude = event.preferred_magnitude() or event.magnitudes[0]

    PreferredTime = origin.time.datetime if origin else None
    PreferredLatitude = origin.latitude if origin else np.nan
    PreferredLongitude = origin.longitude if origin else np.nan
    PreferredDepth = origin.depth / 1000.0 if origin else np.nan  # convert to km
    PreferredMagnitudeValue = magnitude.mag if magnitude else np.nan
    PreferredMagnitudeType = magnitude.magnitude_type if magnitude else ""
    PublicId = event.resource_id.id.split('/')[-1] if event.resource_id else ""

    return {
        "PreferredTime": str(PreferredTime) if PreferredTime is not None else "",
        "PreferredLatitude": float(PreferredLatitude),
        "PreferredLongitude": float(PreferredLongitude),
        "PreferredDepth": float(PreferredDepth),
        "PreferredMagnitudeValue": float(PreferredMagnitudeValue),
        "PreferredMagnitudeType": PreferredMagnitudeType,
        "PublicId": PublicId
    }

if __name__ == '__main__':
    # parsing arguments
    args = parse_args()

    query_and_save(
        starttime=args.starttime,
        endtime=args.endtime,
        minlatitude=args.minlatitude,
        maxlatitude=args.maxlatitude,
        minlongitude=args.minlongitude,
        maxlongitude=args.maxlongitude,
        latitude=args.latitude,
        longitude=args.longitude,
        minradius=args.minradius,
        maxradius=args.maxradius,
        mindepth=args.mindepth,
        maxdepth=args.maxdepth,
        minmagnitude=args.minmagnitude,
        maxmagnitude=args.maxmagnitude,
        eventid=args.eventid,
        output=args.output
    )