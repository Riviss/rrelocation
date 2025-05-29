#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utility script to build GrowClust input files from cross-correlation results.
Edit path constants below before running.

Script to create GrowClust-style files for each cluster:
 - evlist.txt   (date/time lat lon depth placeholders eventID)
 - stlist.txt   (sta lat lon elev)
 - xcordata.txt (combined cross-corr from xcorr_cluster_<cid>_*.txt)

Expects:
 1) A CSV that associates each master_id => (lat, lon, depth, datetime, cluster_id).
 2) A station CSV with columns that include at least the original station code ('sta')
    and a 'unique_id' field that you want to use as the station name.
 3) The xcorr_output/cluster_<cid>/xcorr_cluster_<cid>_*.txt
    files produced by your cross-correlation script with lines like:
      WORKING_cat_A_5367 WORKING_cat_A_1234 BCH1A -0.100 0.965 S
    and with possible header lines starting with '#' that are skipped.

This version maps each master_id to a unique integer, maps station codes to unique_id,
and in the cc (xcorr) files skips self-correlation (i.e. event id 1 = event id 1).
"""

import os
import sys
import pandas as pd
import numpy as np
from datetime import datetime
from obspy import UTCDateTime
import glob

# Centralised path configuration
from config import (
    CLUSTER_ROOT,
    EVENTS_CSV,
    STATIONS_CSV,
)


###############################################################################
# CONFIG
###############################################################################

# Process only clusters >= this value if desired (or leave as None).
MIN_CLUSTER = None  # e.g. 37

# Minimum cross-correlation coefficient to include (0.0 to 1.0)
XCORR_THRESHOLD = 0.6  # Adjust as needed

# Option to square the cross-correlation value before processing (True/False).
# Note: if squared, the threshold applies to the squared value.
SQUARE_CC = False

###############################################################################
# LOAD EVENT / STATION INFO
###############################################################################

def load_event_info(csv_path):
    """
    Loads event data.
    CSV must have columns: master_id, lat, lon, depth, datetime, cluster_id.
    """
    df = pd.read_csv(csv_path)
    df['datetime'] = pd.to_datetime(df['datetime'], errors='coerce')
    return df

def load_station_info(csv_path):
    """
    Loads station info.
    Renames columns if necessary and drops duplicate stations.
    Expected columns after renaming include:
       'sta' (original station code),
       'unique_id' (the desired station name),
       'lat', 'lon', and 'elev'.
    """
    df = pd.read_csv(os.path.expanduser(csv_path))
    df = df.rename({'station': 'sta', 'channel': 'chan', 'network': 'net',
                    'latitude': 'lat', 'longitude': 'lon', 'elevation': 'elev'},
                   axis=1)
    df = df.drop_duplicates(subset='sta', keep='first')
    return df

###############################################################################
# MAPPINGS
###############################################################################

def create_masterid_map(df_events):
    """
    Creates a mapping from master_id (string) to a unique integer.
    """
    unique_ids = df_events['master_id'].unique()
    mapping = {mid: new_id for new_id, mid in enumerate(unique_ids, start=1)}
    return mapping

def create_station_mapping(df_stations):
    """
    Creates a mapping from the original station code ('sta') to the unique station name
    (using the 'unique_id' column, if present). If the unique_id column is missing,
    we simply map to the original value.
    """
    mapping = {}
    if 'unique_id' in df_stations.columns:
        for _, row in df_stations.iterrows():
            mapping[row['sta']] = row['unique_id']
    else:
        for _, row in df_stations.iterrows():
            mapping[row['sta']] = row['sta']
    return mapping

###############################################################################
# 1) CREATE evlist.txt
###############################################################################

def write_evlist(df_events, cluster_id, outdir, master_mapping):
    """
    Writes evlist.txt for a cluster.
    Format example:
      2012 10 13  05 53  3.920  39.66333 -119.68800   7.500    0.01  0.000  0.000  0.000       <unique_event_id>
    """
    outpath = os.path.join(outdir, "evlist.txt")
    skipped_events = 0
    with open(outpath, 'w') as f:
        for _, row in df_events.iterrows():
            master_id_str = row['master_id']
            if master_id_str not in master_mapping:
                print(f"Warning: master_id {master_id_str} not in mapping. Skipping event.")
                skipped_events += 1
                continue
            unique_id = master_mapping[master_id_str]
            lat  = row['lat']
            lon  = row['lon']
            dep  = row['depth']/1000  # m to km
            dt   = row['datetime']
            
            if pd.isna(dt):
                skipped_events += 1
                continue
            
            try:
                year = int(dt.year)
                mon  = int(dt.month)
                day  = int(dt.day)
                hr   = int(dt.hour)
                mn   = int(dt.minute)
                sc   = dt.second + dt.microsecond * 1e-6
            except AttributeError:
                skipped_events += 1
                continue

            # Placeholders as before
            col10, col11, col12, col13 = 0.01, 0.000, 0.000, 0.000

            line = (f"{year:4d} {mon:02d} {day:02d} "
                    f"{hr:02d} {mn:02d} {sc:06.3f} "
                    f"{lat:9.5f} {lon:10.5f} {dep:7.3f} "
                    f"{col10:6.2f} {col11:5.3f} {col12:5.3f} {col13:5.3f} "
                    f"{unique_id:8d}\n")
            f.write(line)
    if skipped_events > 0:
        print(f"  Skipped {skipped_events} events due to errors.")
    print(f"  Wrote evlist.txt => {outpath}")

###############################################################################
# 2) CREATE stlist.txt
###############################################################################

def write_stlist(df_stations, cluster_id, outdir):
    """
    Writes stlist.txt.
    Uses the 'unique_id' from the station CSV if available.
    Format:
      sta lat lon elev(m)
    """
    outpath = os.path.join(outdir, "stlist.txt")
    with open(outpath, 'w') as f:
        # f.write("sta lat lon elev(m)\n")
        for _, row in df_stations.iterrows():
            # Use unique_id if available, else fallback to original 'sta'
            station_name = row['unique_id'] if ('unique_id' in row.index and pd.notnull(row['unique_id'])) else row['sta']
            lat = row['lat']
            lon = row['lon']
            elv = row['elev']
            line = f"{station_name:<5s} {lat:8.4f} {lon:9.4f} {elv:6.1f}\n"
            f.write(line)
    print(f"  Wrote stlist.txt => {outpath}")

###############################################################################
# 3) CREATE xcordata.txt
###############################################################################

def write_xcordata(cluster_id, outdir, master_mapping, station_mapping):
    """
    Reads cross-correlation files with lines like:
      WORKING_cat_A_5367 WORKING_cat_A_1234 BCH1A -0.100 0.965 S
    Produces an aggregated xcordata.txt that:
      - Replaces master_id strings with their unique integer mapping,
      - Remaps the station code to the unique_id (if available), and
      - Skips self-correlations where both event IDs are the same.
    If SQUARE_CC is set to True, the cc value is squared before processing.
    """
    cluster_dir = outdir
    pattern = os.path.join(cluster_dir, f"xcorr_cluster_{cluster_id}_*.txt")
    xcorr_files = glob.glob(pattern)

    from collections import defaultdict
    pair_dict = defaultdict(list)

    for xf in xcorr_files:
        with open(xf, 'r') as f:
            lines = f.read().splitlines()
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Skip header lines (e.g., "# Cross-corr ...")
                continue

            parts = line.split()
            if len(parts) < 6:
                continue

            eventA_str, eventB_str, sta, tdif_str, rxco_str, ph = parts[:6]

            # Check if both events exist in the mapping
            if eventA_str not in master_mapping or eventB_str not in master_mapping:
                print(f"Warning: Event(s) {eventA_str} or {eventB_str} not found in mapping. Skipping line.")
                continue

            idA = master_mapping[eventA_str]
            idB = master_mapping[eventB_str]

            # Skip self-correlation (event compared with itself)
            if idA == idB:
                continue

            try:
                tdif = float(tdif_str)
                rxco = float(rxco_str)
            except ValueError:
                continue

            # Optionally square the cc value
            if SQUARE_CC:
                rxco = rxco * rxco

            # Remap the station name using station_mapping (if available)
            mapped_sta = station_mapping.get(sta, sta)

            # Check threshold using squared value if SQUARE_CC, otherwise use absolute value.
            if SQUARE_CC:
                if rxco >= XCORR_THRESHOLD:
                    pair_key = (idA, idB)
                    pair_dict[pair_key].append((mapped_sta, tdif, rxco, ph))
            else:
                if abs(rxco) >= XCORR_THRESHOLD:
                    pair_key = (idA, idB)
                    pair_dict[pair_key].append((mapped_sta, tdif, rxco, ph))

    outpath = os.path.join(cluster_dir, "xcordata.txt")
    with open(outpath, 'w') as f:
        for (idA, idB) in sorted(pair_dict.keys()):
            f.write(f"#   {idA:6d}   {idB:6d}  0.000\n")
            for entry in pair_dict[(idA, idB)]:
                mapped_sta, tdif, rxco, ph = entry
                f.write(f"    {mapped_sta:<5s} {tdif:8.5f} {rxco:5.4f} {ph}\n")
    print(f"  Wrote xcordata.txt => {outpath}")

###############################################################################
# MAIN
###############################################################################

def main():
    df_events = load_event_info(EVENTS_CSV)
    df_stations = load_station_info(STATIONS_CSV)
    
    # Create mappings for events and stations
    master_mapping = create_masterid_map(df_events)
    station_mapping = create_station_mapping(df_stations)
    print(f"Generated mapping for {len(master_mapping)} master ids.")
    
    cluster_dirs = glob.glob(os.path.join(CLUSTER_ROOT, "cluster_*"))
    print(df_stations)

    for cdir in cluster_dirs:
        basename = os.path.basename(cdir)
        if not basename.startswith("cluster_"):
            continue
        try:
            cid = int(basename.split("_")[-1])
        except ValueError:
            continue

        if MIN_CLUSTER is not None and cid < MIN_CLUSTER:
            continue

        print(f"\n=== Building GrowClust files for {cdir} (cluster {cid}) ===")
        dfc = df_events[df_events['cluster_id'] == cid].copy()
        if dfc.empty:
            print(f"  No events found for cluster {cid}")
            continue

        write_evlist(dfc, cid, cdir, master_mapping)
        write_stlist(df_stations, cid, cdir)
        write_xcordata(cid, cdir, master_mapping, station_mapping)

if __name__ == "__main__":
    main()
