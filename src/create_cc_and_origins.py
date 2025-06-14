#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script prepares cross-correlation files for GrowClust.
Edit the paths in the CONFIG section before running.


DBSCAN + chunked cross-correlation (single-pool approach), reading from CSV files and
waveforms from /home/pgcseiscomp/antelope/wfs/YYYY/MM/DD/YYYYMMDD.*.sta..chan.mseed

Steps:
 1) Load 'origins.csv' & 'arrivals.csv' from your provided paths
 2) DBSCAN on origins => cluster_id
 3) Gather cross-correlation tasks from ALL clusters >= START_CLUSTER, grouped by (network, sta, chan, phase)
 4) Process all tasks with one multiprocessing pool. Each task:
    - chunk-splits picks so we never load huge waveforms at once
    - loads waveforms (bandpass+trim) with a wildcarded network code
    - cross-correlates pairs
    - batch-writes lines to a .txt file so we avoid large returns.

We skip cluster IDs below START_CLUSTER. Log warnings/errors for missing waveforms, etc.
"""

import os
import glob
import numpy as np
import pandas as pd
import warnings
import multiprocessing
import argparse
import logging
from multiprocessing import Pool
from datetime import datetime

from obspy import read, UTCDateTime
from obspy.signal.cross_correlation import correlate
from sklearn.cluster import DBSCAN

# Centralised path configuration
# Import configuration from the src package. We use an absolute import so that
# this script can be executed directly as ``python src/create_cc_and_origins.py``
# without requiring ``-m`` module syntax.
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
    WAVEFORM_ROOT,
    XCORR_OUTDIR,
    get_db_engine,
)

###############################################################################
# Logging
###############################################################################

def setup_logging(level: str = "INFO") -> None:
    """Configure root logger."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

warnings.filterwarnings("ignore", category=DeprecationWarning)

###############################################################################
# OUTPUT CONFIG
###############################################################################
os.makedirs(XCORR_OUTDIR, exist_ok=True)

###############################################################################
# DBSCAN SETTINGS
###############################################################################
KM_RADIUS       = 1.0   # 'eps' in XY
TIME_WINDOW_DAYS= 7.0   # scale 7 days => 1 km in t-dimension
MIN_SAMPLES     = 10    # minimum # of events per cluster
MAX_CLUSTER_EVENTS = 1000  # recursively split clusters larger than this

###############################################################################
# CROSS-CORRELATION SETTINGS
###############################################################################
PRE_TIME       = 0.1      # seconds before pick
POST_TIME      = 0.4      # seconds after pick
BUFFER_SEC     = 1.0      # initial trim buffer
EDGE_DROP_SEC  = 0.2      # drop edges after filtering
CHUNK_SIZE     = 100      # picks per chunk => avoid memory blow-up
MAX_SHIFT_SEC  = 0.4      # ±0.4 s shift
CC_THRESHOLD   = 0.6      # |ρ| ≥ threshold
ABSOLUTE_CORR  = True     # if True => require abs(ccmax) >= CC_THRESHOLD
BANDPASS_LO    = 2.0
BANDPASS_HI    = 10.0

# We'll skip any cluster with ID < START_CLUSTER
START_CLUSTER  = 0

###############################################################################
# LOAD + DBSCAN
###############################################################################
def approximate_xy(lat, lon, lat0=None, lon0=None):
    EARTH_RADIUS_KM = 6371.0
    DEG2RAD         = np.pi/180.0
    if lat0 is None:
        lat0 = np.mean(lat)
    if lon0 is None:
        lon0 = np.mean(lon)
    x = EARTH_RADIUS_KM * np.cos(lat0*DEG2RAD)*(lon - lon0)*DEG2RAD
    y = EARTH_RADIUS_KM*(lat - lat0)*DEG2RAD
    return x, y

def load_and_cluster_dbscan():
    """
    1) Load origin and arrival information
       either from CSV files or from the configured MySQL database.
    2) DBSCAN on origins => cluster_id
    3) Return df_origin, df_arrival (not merged yet)
    """

    if USE_DATABASE:
        engine = get_db_engine()
        df_origin = pd.read_sql(f"SELECT * FROM {ORIGINS_TABLE}", engine)
        df_arrival = pd.read_sql(f"SELECT * FROM {ARRIVALS_TABLE}", engine)
    else:
        df_origin = pd.read_csv(ORIGINS_PATH)
        df_arrival = pd.read_csv(ARRIVALS_PATH)

    # rename origin columns to standard
    df_origin.rename(columns={
        'event_datetime': 'datetime',
        'latitude': 'lat',
        'longitude': 'lon',
        'depth_m': 'depth'
    }, inplace=True, errors='ignore')
    df_origin['datetime'] = pd.to_datetime(df_origin['datetime'], errors='coerce')
    df_origin.dropna(subset=['master_id','datetime','lat','lon'], inplace=True)

    # rename arrival columns
    if 'datetime' in df_arrival.columns and 'arrival_time' not in df_arrival.columns:
        df_arrival.rename(columns={'datetime':'arrival_time'}, inplace=True)
    df_arrival['arrival_time'] = pd.to_datetime(df_arrival['arrival_time'], errors='coerce')
    df_arrival.dropna(subset=['master_id','arrival_time','sta','phase'], inplace=True)

    # DBSCAN in (x,y,t) with recursive splitting for large clusters
    lat_arr = df_origin['lat'].values
    lon_arr = df_origin['lon'].values
    x_arr, y_arr = approximate_xy(lat_arr, lon_arr)

    dt_min = df_origin['datetime'].min()
    dt_days = (df_origin['datetime'] - dt_min).dt.total_seconds() / 86400.0
    scale_time = KM_RADIUS / TIME_WINDOW_DAYS
    t_scaled = dt_days * scale_time

    coords = np.column_stack([x_arr, y_arr, t_scaled])

    def recursive_dbscan(idxs):
        nonlocal next_id
        subcoords = coords[idxs]
        labels = DBSCAN(eps=KM_RADIUS, min_samples=MIN_SAMPLES,
                        metric='chebyshev', n_jobs=-1).fit_predict(subcoords)
        for lbl in np.unique(labels):
            mask = labels == lbl
            sel = idxs[mask]
            if lbl == -1:
                final_labels[sel] = -1
            elif len(sel) <= MAX_CLUSTER_EVENTS:
                final_labels[sel] = next_id
                next_id += 1
            else:
                stack.append(sel)

    n = len(coords)
    final_labels = np.full(n, -1, dtype=int)
    next_id = 0
    stack = [np.arange(n)]
    while stack:
        recursive_dbscan(stack.pop())

    df_origin['cluster_id'] = final_labels.astype(int)
    return df_origin, df_arrival

###############################################################################
# LOAD WAVEFORM FOR SINGLE PICK
###############################################################################
def load_trace_for_pick(row):
    """
    Build path using wildcard for network code:
      /home/pgcseiscomp/antelope/wfs/YYYY/MM/DD/YYYYMMDD.*.STA..CHAN.mseed
    Then read, bandpass, final trim => return single Trace or None.
    """
    pick_t = row['pick_time']  # an ObsPy UTCDateTime
    yyyy   = pick_t.year
    mm     = pick_t.month
    dd     = pick_t.day
    ymd    = pick_t.strftime("%Y%m%d")

    # We'll ignore 'network' or 'snet' in the filename => wildcard
    sta = row['sta']
    cha = row['chan']

    subdir = os.path.join(WAVEFORM_ROOT, f"{yyyy:04d}", f"{mm:02d}", f"{dd:02d}")
    pattern = f"{ymd}.*.{sta}..{cha}.mseed"
    path_pattern = os.path.join(subdir, pattern)

    matches = glob.glob(path_pattern)
    if not matches:
        logging.warning("no file match => %s", path_pattern)
        return None

    f = matches[0]
    try:
        st = read(f)
        # 1) wide trim with buffer
        st.trim(pick_t - PRE_TIME - BUFFER_SEC, pick_t + POST_TIME + BUFFER_SEC)
        # 2) zero-phase bandpass and taper
        st.filter(
            "bandpass",
            freqmin=BANDPASS_LO,
            freqmax=BANDPASS_HI,
            corners=4,
            zerophase=True,
        )
        st.taper(max_percentage=0.05, type="cosine")
        # 3) demean and detrend
        st.detrend("demean")
        st.detrend("linear")
        # 4) drop filter edge before final window
        st.trim(
            pick_t - PRE_TIME - EDGE_DROP_SEC,
            pick_t + POST_TIME + EDGE_DROP_SEC,
        )
        # 5) final precise trim
        st.trim(pick_t - PRE_TIME, pick_t + POST_TIME)

        if len(st) == 0:
            logging.warning("no data after trim => %s, pick@%s", f, pick_t)
            return None
        tr = st[0]
        tr.stats.master_id = row['master_id']
        return tr
    except Exception as e:
        logging.error("loading %s failed: %s", f, e)
        return None

###############################################################################
# NORMALIZED CROSS-CORRELATION HELPER
###############################################################################
def ncc_max(dataA: np.ndarray, dataB: np.ndarray, max_samps: int) -> tuple[float, int]:
    """Return the maximum normalized correlation and lag in samples."""
    cc = correlate(dataA, dataB, max_samps)
    norm = np.linalg.norm(dataA) * np.linalg.norm(dataB)
    if norm == 0:
        return 0.0, 0
    cc /= norm
    imax = np.argmax(np.abs(cc))
    return float(cc[imax]), imax - max_samps

###############################################################################
# CROSS-CORRELATION (CHUNKED) WORKER
###############################################################################
def crosscorr_group(args):
    """
    Worker function for a single (network, sta, chan, phase) group in a cluster.
    We chunk-split picks => load waveforms => cross-corr => write .txt
    """
    cid    = str(int(args['cluster_id']))
    net    = args['network']
    sta    = args['sta']
    cha    = args['chan']
    phase  = args['phase']
    df_g   = args['df_group']
    cdir   = args['cluster_dir']

    outname = f"xcorr_cluster_{cid}_{net}_{sta}_{cha}_{phase}.txt"
    outpath = os.path.join(cdir, outname)

    logging.info(
        "[cluster %s] crosscorr_group net=%s sta=%s chan=%s phase=%s picks=%d",
        cid,
        net,
        sta,
        cha,
        phase,
        len(df_g),
    )

    with open(outpath, 'w') as f_out:
        hdr = f"# Cross-corr cluster={cid}, network={net}, sta={sta}, chan={cha}, phase={phase}, picks={len(df_g)}"
        f_out.write(hdr + "\n")

        if len(df_g) < 2:
            f_out.write("# Not enough picks to form pairs\n")
            return f"Cluster {cid}:{net}.{sta}.{cha}.{phase} => <2 picks => done."

        # Sort by pick_time if desired
        df_g = df_g.sort_values('pick_time').reset_index(drop=True)
        n = len(df_g)
        n_chunks = int(np.ceil(n / CHUNK_SIZE))

        lines_buffer = []
        BATCH_SIZE   = 500

        def flush_buffer():
            if lines_buffer:
                f_out.write("\n".join(lines_buffer) + "\n")
                lines_buffer.clear()
                f_out.flush()

        def load_chunk(df_subset):
            waveforms = []
            for _, rowp in df_subset.iterrows():
                tr = load_trace_for_pick(rowp)
                if tr is not None:
                    waveforms.append(tr)
            return waveforms

        for i in range(n_chunks):
            i0 = i * CHUNK_SIZE
            i1 = min((i+1)*CHUNK_SIZE, n)
            df_i = df_g.iloc[i0:i1]
            st_i = load_chunk(df_i)

            # Self correlation
            if len(st_i) > 1:
                for ii in range(len(st_i)-1):
                    trA = st_i[ii]
                    midA= trA.stats.master_id
                    srA = trA.stats.sampling_rate
                    dataA= trA.data
                    for jj in range(ii+1, len(st_i)):
                        trB = st_i[jj]
                        midB= trB.stats.master_id
                        srB = trB.stats.sampling_rate
                        dataB= trB.data

                        # unify sample rate
                        if abs(srA - srB) > 0.3:
                            trB = trB.copy()
                            trB.resample(srA)
                            dataB= trB.data
                            sr= srA
                        else:
                            sr= srA

                        max_samps = int(round(MAX_SHIFT_SEC * sr))
                        ccmax, lag_samp = ncc_max(dataA, dataB, max_samps)
                        if abs(ccmax) < CC_THRESHOLD:
                            continue
                        lag_sec  = lag_samp / sr
                        line = f"{midA} {midB} {sta} {lag_sec:.3f} {ccmax:.3f} {phase}"
                        lines_buffer.append(line)
                        if len(lines_buffer) >= BATCH_SIZE:
                            flush_buffer()

            # cross-chunk st_i vs st_j
            for j in range(i+1, n_chunks):
                j0 = j * CHUNK_SIZE
                j1 = min((j+1)*CHUNK_SIZE, n)
                df_j = df_g.iloc[j0:j1]
                st_j = load_chunk(df_j)
                if len(st_j)==0 or len(st_i)==0:
                    st_j.clear()
                    continue

                for trA in st_i:
                    midA= trA.stats.master_id
                    srA = trA.stats.sampling_rate
                    dataA= trA.data
                    for trB in st_j:
                        midB= trB.stats.master_id
                        srB = trB.stats.sampling_rate
                        dataB= trB.data

                        if abs(srA - srB) > 0.3:
                            trB= trB.copy()
                            trB.resample(srA)
                            dataB= trB.data
                            sr= srA
                        else:
                            sr= srA

                        max_samps = int(round(MAX_SHIFT_SEC * sr))
                        ccmax, lag_samp = ncc_max(dataA, dataB, max_samps)
                        if abs(ccmax) < CC_THRESHOLD:
                            continue
                        lag_sec  = lag_samp / sr
                        line = f"{midA} {midB} {sta} {lag_sec:.3f} {ccmax:.3f} {phase}"
                        lines_buffer.append(line)
                        if len(lines_buffer) >= BATCH_SIZE:
                            flush_buffer()

                st_j.clear()

            st_i.clear()

        flush_buffer()

    return f"Cluster {cid}:{net}.{sta}.{cha}.{phase} => wrote {outname}"

###############################################################################
# MAIN (SINGLE POOL FOR ALL TASKS)
###############################################################################
def run_script(start_cluster: int = START_CLUSTER,
               cc_threshold: float = CC_THRESHOLD,
               processes: int | None = None) -> None:
    """Main entry point for the cross-correlation workflow."""
    global START_CLUSTER, CC_THRESHOLD

    START_CLUSTER = start_cluster
    CC_THRESHOLD = cc_threshold

    cpu_count = processes or multiprocessing.cpu_count()
    logging.info("Detected %d CPU cores. Using %d processes", cpu_count, cpu_count)

    # 1) Load & DBSCAN
    df_origin, df_arrival = load_and_cluster_dbscan()
    out_csv = os.path.join(XCORR_OUTDIR, "clustered_origins.csv")
    df_origin.to_csv(out_csv, index=False)
    logging.info("DBSCAN done. Wrote => %s", out_csv)

    # 2) Merge cluster_id => arrivals, skip noise
    df_arrival = pd.merge(df_arrival, df_origin[['master_id','cluster_id']],
                          how='left', on='master_id')
    df_arrival = df_arrival[df_arrival['cluster_id']!=-1].copy()

    # Convert arrival_time => pick_time
    df_arrival['pick_time'] = df_arrival['arrival_time'].apply(lambda x: UTCDateTime(x))

    # (Optional) load stations => 'network' column
    if os.path.exists(STATIONS_PATH):
        df_st = pd.read_csv(STATIONS_PATH)[['unique_id','station','network']].drop_duplicates()
        df_st.rename(columns={'station':'sta'}, inplace=True)
        df_arrival = pd.merge(df_arrival, df_st, on='sta', how='left')
        df_arrival['network'] = df_arrival['network'].fillna('NA')
    else:
        df_arrival['network'] = 'NA'

    # 3) Gather tasks for ALL clusters
    tasks = []
    clusters = sorted(df_arrival['cluster_id'].unique().astype(int))
    logging.info(
        "Found %d clusters (excl noise). Clusters: %s",
        len(clusters),
        clusters,
    )

    for cid in clusters:
        if cid < START_CLUSTER:
            continue
        dfc = df_arrival[df_arrival['cluster_id']==cid]
        n_events = df_origin[df_origin['cluster_id']==cid].shape[0]
        if n_events < 2:
            continue

        # Each cluster => subdir
        cdir = os.path.join(XCORR_OUTDIR, f"cluster_{cid}")
        os.makedirs(cdir, exist_ok=True)

        # group by (network, sta, chan, phase)
        group_cols = ['network','sta','chan','phase']
        for (net, st, ch, ph), df_g in dfc.groupby(group_cols):
            if len(df_g)<2:
                continue
            tasks.append({
                'cluster_id': cid,
                'network':   net,
                'sta':       st,
                'chan':      ch,
                'phase':     ph,
                'df_group':  df_g,
                'cluster_dir': cdir
            })

    logging.info(
        "Prepared %d cross-correlation tasks across all clusters.", len(tasks)
    )

    if tasks:
        # 4) Single pool => keep all cores busy
        with Pool(processes=cpu_count) as pool:
            results = pool.map(crosscorr_group, tasks)

        for msg in results:
            logging.info("Worker => %s", msg)
    else:
        logging.warning("No cross-correlation tasks found!")

    logging.info("All done. Check xcorr_output/cluster_* directories.")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Prepare cross-correlation files for GrowClust"
    )
    parser.add_argument(
        "--start-cluster",
        type=int,
        default=START_CLUSTER,
        help="Skip clusters with ID below this value",
    )
    parser.add_argument(
        "--cc-threshold",
        type=float,
        default=CC_THRESHOLD,
        help="Minimum correlation coefficient to keep a pair",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Number of worker processes to use",
    )
    parser.add_argument(
        "--log",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR)",
    )

    args = parser.parse_args()
    setup_logging(args.log)
    run_script(
        start_cluster=args.start_cluster,
        cc_threshold=args.cc_threshold,
        processes=args.processes,
    )
