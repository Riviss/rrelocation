#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chunk-based cross-correlation with minimal memory usage:
Edit the CONFIG section below before running this script.

 - We only ever load two chunks of waveforms (chunk i and chunk j)
   into memory at a time.
 - Immediately write cross-correlation results to CSV to avoid
   storing large arrays in memory.
"""

import os
import glob
import gc
import numpy as np
import pandas as pd

from obspy import read, Stream
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal import cross_correlation as xcorr

from multiprocessing import Pool
from sklearn.cluster import DBSCAN

# Centralised path configuration
from config import (
    USE_DATABASE,
    DB_HOST,
    DB_USER,
    DB_PASSWORD,
    DB_NAME,
    ORIGINS_TABLE,
    ARRIVALS_TABLE,
    ORIGINS_PATH,
    ARRIVALS_PATH,
    STATIONS_PATH,
    WAVEFORM_ROOT as WF_PATH,
    XCORR_OUTDIR,
    get_db_engine,
)

EARTH_RADIUS_KM = 6371.0
DEG2RAD         = np.pi / 180.0

TIME_WINDOW_DAYS   = 7
DIST_KM            = 2
DBSCAN_MIN_SAMPLES = 10
MAX_CLUSTER_EVENTS = 1000

PRETRIG_DEFAULT  = 0.1
POSTTRIG_DEFAULT = 0.4
LEAD_TIME_DEFAULT = 3
POST_TIME_DEFAULT = 3

NUM_CORES  = 8

# If memory is tight, reduce this. Each chunk will hold CHUNK_SIZE waveforms max.
CHUNK_SIZE = 200

# Preprocessing parameters
BANDPASS_LO = 2.0  # Hz
BANDPASS_HI = 10.0  # Hz

###############################################################################
# DBSCAN & Data-Loading Helpers
###############################################################################

def approximate_xy(lat, lon, lat0=None, lon0=None):
    """Convert lat/lon to approximate local (X,Y) coordinates in km."""
    if lat0 is None:
        lat0 = np.mean(lat)
    if lon0 is None:
        lon0 = np.mean(lon)
    x = EARTH_RADIUS_KM * np.cos(lat0 * DEG2RAD) * (lon - lon0) * DEG2RAD
    y = EARTH_RADIUS_KM * (lat - lat0) * DEG2RAD
    return x, y

def cluster_events_dbscan(df_origins, time_window_days=7, dist_km=2, min_samples=1):
    """
    Cluster events in a 3D space (X, Y, T_scaled) using DBSCAN (Chebyshev metric)
    and recursively split clusters larger than ``MAX_CLUSTER_EVENTS``.
    """
    df = df_origins.copy().sort_values(by='datetime').reset_index(drop=True)
    lat_arr = df['lat'].values
    lon_arr = df['lon'].values
    x_arr, y_arr = approximate_xy(lat_arr, lon_arr)
    dt_min = df['datetime'].min()
    dt_days = (df['datetime'] - dt_min).dt.total_seconds() / 86400.0

    time_scale = dist_km / time_window_days
    t_scaled = dt_days * time_scale

    coords = np.column_stack([x_arr, y_arr, t_scaled])

    def recursive_dbscan(idxs):
        nonlocal next_id
        subcoords = coords[idxs]
        labels = DBSCAN(eps=dist_km, min_samples=min_samples,
                        metric='chebyshev', n_jobs=-1).fit_predict(subcoords)
        for lbl in np.unique(labels):
            subset = idxs[labels == lbl]
            if lbl == -1:
                final_labels[subset] = -1
            elif len(subset) <= MAX_CLUSTER_EVENTS:
                final_labels[subset] = next_id
                next_id += 1
            else:
                stack.append(subset)

    n = len(coords)
    final_labels = np.full(n, -1, dtype=int)
    next_id = 0
    stack = [np.arange(n)]
    while stack:
        recursive_dbscan(stack.pop())

    df['cluster'] = final_labels.astype(int)
    return df

def load_and_cluster_dbscan():
    """
    Load origin and arrival information from either CSV files or the
    configured MySQL database and perform DBSCAN clustering.
    Returns (origin_df, arrival_df).
    """
    if USE_DATABASE:
        engine = get_db_engine()
        origin = pd.read_sql(f"SELECT * FROM {ORIGINS_TABLE}", engine)
        arrival = pd.read_sql(f"SELECT * FROM {ARRIVALS_TABLE}", engine)
    else:
        origin = pd.read_csv(ORIGINS_PATH)
        arrival = pd.read_csv(ARRIVALS_PATH)
    
    # Rename columns to standard
    origin.rename(columns={
        'event_datetime': 'datetime',
        'latitude': 'lat',
        'longitude': 'lon',
        'depth_m': 'depth'
    }, inplace=True)
    
    origin['datetime'] = pd.to_datetime(origin['datetime'], errors='coerce')
    origin.dropna(subset=['datetime', 'lat', 'lon'], inplace=True)
    
    arrival['datetime'] = pd.to_datetime(arrival['datetime'], errors='coerce')
    arrival.dropna(subset=['datetime'], inplace=True)
    
    # DBSCAN
    origin = cluster_events_dbscan(origin,
                                   time_window_days=TIME_WINDOW_DAYS,
                                   dist_km=DIST_KM,
                                   min_samples=DBSCAN_MIN_SAMPLES)
    
    return origin, arrival

def refine_data(origin, arrival,
                pretrig=0.3, posttrig=0.5, lead_time=3, post_time=3,
                datetime_col='datetime', station_col='sta'):
    """
    Merge station info, build waveform paths, add start/end times to arrival, etc.
    """
    stations = pd.read_csv(STATIONS_PATH)[['unique_id', 'station']]
    
    arrival = arrival.copy()
    arrival = pd.merge(arrival, stations, how='left', left_on='sta', right_on='station')
    arrival['snet'] = arrival['unique_id']
    
    origin['ortime'] = origin[datetime_col].apply(lambda x: UTCDateTime(x))
    arrival['artime'] = arrival[datetime_col].apply(lambda x: UTCDateTime(x))
    
    arrival['pretrig'] = pretrig
    arrival['posttrig'] = posttrig
    
    arrival['artime_start'] = arrival['artime'] - lead_time - pretrig
    arrival['artime_end']   = arrival['artime'] + post_time  + posttrig

    arrival['path'] = (
        WF_PATH +
        arrival[datetime_col].dt.strftime('%Y/%m/%d/') +
        arrival[datetime_col].dt.strftime('%Y%m%d') +
        ".*." + arrival[station_col].fillna('NA') +
        ".." + arrival['chan'].fillna('NA') +
        ".mseed"
    )
    
    return origin, arrival

###############################################################################
# Waveform & Xcorr Helpers
###############################################################################

def load_waveforms_for_arrivals(df):
    """
    Given a chunk of arrival rows, read & trim each waveform from disk,
    and return as an ObsPy Stream. Each row => 1 Trace if found.
    """
    st = Stream()
    for _, row in df.iterrows():
        pattern = row["path"]
        matches = glob.glob(pattern)
        if not matches:
            continue
        f = matches[0]
        try:
            wf = read(f)
            # initial window with buffer
            wf.trim(row["artime_start"], row["artime_end"])
            wf.filter(
                "bandpass",
                freqmin=BANDPASS_LO,
                freqmax=BANDPASS_HI,
                corners=4,
                zerophase=True,
            )
            wf.detrend("demean")
            wf.detrend("linear")
            # final precise window after filtering
            wf.trim(
                row["artime"] - row["pretrig"],
                row["artime"] + row["posttrig"],
            )
            if len(wf) > 0:
                tr = wf[0]
                tr.stats.station = str(row["sta"])
                tr.stats.channel = str(row["chan"])
                tr.stats.phase = str(row["phase"])
                tr.stats.artime = row["artime"]
                tr.stats.orid = row["master_id"]
                st.append(tr)
        except Exception as e:
            print(f"Error reading {f}: {e}")
    return st

def xcorr_pairwise(stA, stB, max_lag=50, self_xcorr=False):
    """
    Cross-correlate all pairs in stA × stB. If self_xcorr=True, stA==stB
    and we only do i<j. Returns list of dict results.
    """
    results = []
    if self_xcorr:
        for i in range(len(stA) - 1):
            for j in range(i + 1, len(stA)):
                trA = stA[i]
                trB = stA[j]
                cc = xcorr.correlate(trA.data, trB.data, max_lag)
                norm = np.linalg.norm(trA.data) * np.linalg.norm(trB.data)
                if norm != 0:
                    cc /= norm
                    cc_max = cc[np.argmax(np.abs(cc))]
                else:
                    cc_max = 0.0
                results.append({
                    "oridA": trA.stats.orid,
                    "oridB": trB.stats.orid,
                    "station": trA.stats.station,
                    "channel": trA.stats.channel,
                    "phase": trA.stats.phase,
                    "xcorr": cc_max,
                })
    else:
        for trA in stA:
            for trB in stB:
                cc = xcorr.correlate(trA.data, trB.data, max_lag)
                norm = np.linalg.norm(trA.data) * np.linalg.norm(trB.data)
                if norm != 0:
                    cc /= norm
                    cc_max = cc[np.argmax(np.abs(cc))]
                else:
                    cc_max = 0.0
                results.append({
                    "oridA": trA.stats.orid,
                    "oridB": trB.stats.orid,
                    "station": trA.stats.station,
                    "channel": trA.stats.channel,
                    "phase": trA.stats.phase,
                    "xcorr": cc_max,
                })
    return results

###############################################################################
# Cross-Correlation Task (chunk-based, 2-layer)
###############################################################################

def crosscorr_station_phase(args):
    """
    For one station-phase-chunk set:
      - We have a subset of arrivals 'df' for this station/chan/phase/cluster.
      - We chunk them, load chunk i as st_i, chunk j as st_j, cross-correlate, write results,
        discard st_i/st_j. This drastically limits memory usage.
    """
    cid = args['cluster_id']
    sta = args['station']
    cha = args['channel']
    pha = args['phase']
    df  = args['arr_subset']
    
    # We'll store results in a single CSV for that station–phase:
    out_file = os.path.join(XCORR_OUTDIR,
                            f"cluster_{cid}_station_{sta}_phase_{pha}.csv")
    with open(out_file, "w") as f:
        f.write("oridA,oridB,station,channel,phase,xcorr\n")

    # Track written pairs to avoid duplicates.  We normalise the event order
    # (min, max) so A,B and B,A are considered the same pair.
    written_pairs = set()
    
    # If fewer than 2 arrivals, skip
    n = len(df)
    if n < 2:
        print(f"Cluster {cid} {sta}.{cha}.{pha}: Only {n} arrivals => skip")
        return
    
    # Build chunk boundaries
    n_chunks = int(np.ceil(n / CHUNK_SIZE))
    
    # Indices for chunk i
    for i in range(n_chunks):
        i0 = i * CHUNK_SIZE
        i1 = min((i + 1) * CHUNK_SIZE, n)
        df_i = df.iloc[i0:i1]
        st_i = load_waveforms_for_arrivals(df_i)

        # Self cross-correlation for chunk i
        if len(st_i) > 1:
            results_i = xcorr_pairwise(st_i, st_i, max_lag=50, self_xcorr=True)
            if results_i:
                filtered = []
                for r in results_i:
                    key = tuple(sorted((r["oridA"], r["oridB"]))) + (
                        r["station"], r["channel"], r["phase"])
                    if key not in written_pairs:
                        written_pairs.add(key)
                        filtered.append(r)
                if filtered:
                    pd.DataFrame(filtered).to_csv(
                        out_file, mode="a", header=False, index=False)

        # Cross-correlation with subsequent chunk j>i
        for j in range(i+1, n_chunks):
            j0 = j * CHUNK_SIZE
            j1 = min((j + 1) * CHUNK_SIZE, n)
            df_j = df.iloc[j0:j1]
            st_j = load_waveforms_for_arrivals(df_j)
            
            if len(st_i) == 0 or len(st_j) == 0:
                del st_j
                gc.collect()
                continue
            
            # cross correlate st_i vs st_j
            results_ij = xcorr_pairwise(st_i, st_j, max_lag=50, self_xcorr=False)
            if results_ij:
                filtered = []
                for r in results_ij:
                    key = tuple(sorted((r["oridA"], r["oridB"]))) + (
                        r["station"], r["channel"], r["phase"])
                    if key not in written_pairs:
                        written_pairs.add(key)
                        filtered.append(r)
                if filtered:
                    pd.DataFrame(filtered).to_csv(
                        out_file, mode="a", header=False, index=False)
            
            # discard st_j
            del st_j
            gc.collect()
        
        # discard st_i
        del st_i
        gc.collect()

    print(f"Finished cross-correlation for cluster {cid}, station {sta}, phase {pha}")

###############################################################################
# Main script
###############################################################################

def run_cc():
    print("=== Loading & clustering data ===")
    origin, arrival = load_and_cluster_dbscan()
    print(f"Loaded {len(origin)} origins, {len(arrival)} arrivals.")
    print("Clusters found:", origin['cluster'].unique())

    # Save cluster results
    out_clust = "clustered_origins.csv"
    origin.to_csv(out_clust, index=False)
    print(f"Clustered origin data saved to {out_clust}")

    # Refine
    origin_ref, arrival_ref = refine_data(origin, arrival,
                                          pretrig=PRETRIG_DEFAULT,
                                          posttrig=POSTTRIG_DEFAULT,
                                          lead_time=LEAD_TIME_DEFAULT,
                                          post_time=POST_TIME_DEFAULT)
    print("Refined data => ", len(origin_ref), "origins,", len(arrival_ref), "arrivals.")

    # Keep only real clusters (ignore cluster == -1)
    valid_clusters = sorted(origin_ref[origin_ref['cluster'] != -1]['cluster'].unique())
    print(f"Valid clusters: {valid_clusters}")

    for cid in valid_clusters:
        cluster_mask = (origin_ref['cluster'] == cid)
        c_events = origin_ref.loc[cluster_mask, 'master_id'].unique()
        if len(c_events) < 2:
            print(f"Skipping cluster {cid} (only {len(c_events)} events).")
            continue

        # Subset arrivals for just these events
        arr_c = arrival_ref[arrival_ref['master_id'].isin(c_events)]
        if arr_c.empty:
            print(f"No arrivals found for cluster {cid}, skip.")
            continue

        # Group by station, channel, phase
        tasks = []
        for (sta, cha, pha), df_sub in arr_c.groupby(['sta', 'chan', 'phase']):
            if len(df_sub) < 2:
                continue
            tasks.append({
                'cluster_id': cid,
                'station': sta,
                'channel': cha,
                'phase': pha,
                'arr_subset': df_sub
            })
        
        print(f"Cluster {cid} => {len(tasks)} station-phase tasks")
        if not tasks:
            continue

        # Parallelize
        with Pool(processes=min(NUM_CORES, len(tasks))) as pool:
            for _ in pool.imap_unordered(crosscorr_station_phase, tasks, chunksize=1):
                pass
        
        del tasks, arr_c
        gc.collect()

    print("=== All cross-correlations complete ===")

if __name__ == "__main__":
    run_cc()
