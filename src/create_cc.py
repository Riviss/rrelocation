#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Robust, memory‑efficient cross‑correlation script for GrowClust / hypocentral relocation.

Features added to the original `create_cc_and_origins.py`:
  • **Two‑stage trimming + taper + demean + detrend** for clean CC windows
  • **Normalized cross‑correlation (NCC)** with absolute‑value threshold
  • **CSV output** (one file per cluster×station×phase)
  • **Dynamic chunk size** and multiprocessing
  • **Extensive CLI flags** and logging
  • Keeps all DBSCAN, database/CSV IO, and task‑assembly logic from the original

Usage example:
    python create_cc_refined.py \
        --bandpass-lo 1.5 --bandpass-hi 12 \
        --pre 0.1 --post 0.4 \
        --buffer 1.0 --edge-drop 0.2 \
        --cc-threshold 0.6 \
        --chunk-size 150 \
        --processes 16
"""

###############################################################################
# Imports & configuration
###############################################################################
import os
import glob
import gc
import argparse
import logging
import multiprocessing as mp
from multiprocessing import Pool
from datetime import datetime

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

from obspy import read, UTCDateTime, Stream
from obspy.signal.cross_correlation import correlate

from config import (
    USE_DATABASE,
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
# Default parameters (over‑ridable from CLI)
###############################################################################
KM_RADIUS_DEF       = 1.0      # km → DBSCAN eps in X,Y
TIME_WINDOW_DAYS_DEF= 7.0      # days → 1 km in time dimension
MIN_SAMPLES_DEF     = 10

PRE_TIME_DEF        = 0.1      # s before pick (final window)
POST_TIME_DEF       = 0.4      # s after pick (final window)
BUFFER_DEF          = 1.0      # s extra padding around initial cut
EDGE_DROP_DEF       = 0.2      # s dropped after filter (before final cut)

BANDPASS_LO_DEF     = 2.0      # Hz
BANDPASS_HI_DEF     = 10.0     # Hz

MAX_SHIFT_SEC_DEF   = 0.4      # ± shift for CC (s)
CC_THRESHOLD_DEF    = 0.6      # |ρ| ≥ this keeps pair
CHUNK_SIZE_DEF      = 100      # arrivals per chunk

START_CLUSTER_DEF   = 0        # skip clusters with ID < this
LOGLEVEL_DEF        = "INFO"

###############################################################################
# Helpers
###############################################################################

def setup_logging(level:str):
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def approximate_xy(lat, lon, lat0=None, lon0=None):
    """Quick lat/lon → local X,Y (km) approximation."""
    R = 6371.0
    d2r = np.pi / 180.0
    lat0 = np.mean(lat) if lat0 is None else lat0
    lon0 = np.mean(lon) if lon0 is None else lon0
    x = R * np.cos(lat0 * d2r) * (lon - lon0) * d2r
    y = R * (lat - lat0) * d2r
    return x, y


def load_and_cluster_dbscan(km_radius, time_win_days, min_samples):
    """Load origins/arrivals (CSV or DB) and run 3‑D DBSCAN."""
    if USE_DATABASE:
        engine = get_db_engine()
        df_o = pd.read_sql(f"SELECT * FROM {ORIGINS_TABLE}", engine)
        df_a = pd.read_sql(f"SELECT * FROM {ARRIVALS_TABLE}", engine)
    else:
        df_o = pd.read_csv(ORIGINS_PATH)
        df_a = pd.read_csv(ARRIVALS_PATH)

    df_o.rename(columns={
        "event_datetime": "datetime",
        "latitude": "lat",
        "longitude": "lon",
        "depth_m": "depth",
    }, inplace=True, errors="ignore")
    df_o["datetime"] = pd.to_datetime(df_o["datetime"], errors="coerce")
    df_o.dropna(subset=["master_id", "datetime", "lat", "lon"], inplace=True)

    if "datetime" in df_a.columns and "arrival_time" not in df_a.columns:
        df_a.rename(columns={"datetime": "arrival_time"}, inplace=True)
    df_a["arrival_time"] = pd.to_datetime(df_a["arrival_time"], errors="coerce")
    df_a.dropna(subset=["master_id", "arrival_time", "sta", "phase"], inplace=True)

    x, y = approximate_xy(df_o["lat"].values, df_o["lon"].values)
    dt0 = df_o["datetime"].min()
    t_days = (df_o["datetime"] - dt0).dt.total_seconds() / 86400.0
    t_scaled = t_days * km_radius / time_win_days

    coords = np.column_stack([x, y, t_scaled])
    db = DBSCAN(eps=km_radius, min_samples=min_samples, metric="chebyshev", n_jobs=-1)
    df_o["cluster_id"] = db.fit_predict(coords).astype(int)

    return df_o, df_a

###############################################################################
# Waveform loading & preprocessing
###############################################################################

def preprocess_trace(fpath:str, pick_t:UTCDateTime, pre:float, post:float,
                     buf:float, edge_drop:float, f_lo:float, f_hi:float,
                     sta:str, cha:str, master_id:int):
    try:
        st = read(fpath)
        # 1⃣ wide cut with buffer
        st.trim(pick_t - pre - buf, pick_t + post + buf)
        # 2⃣ band‑pass (zero‑phase 4‑pole) + taper
        st.filter("bandpass", freqmin=f_lo, freqmax=f_hi, corners=4, zerophase=True)
        st.taper(max_percentage=0.05, type="cosine")
        # 3⃣ demean & detrend
        st.detrend("demean")
        st.detrend("linear")
        # 4⃣ drop filter edge
        st.trim(pick_t - pre - edge_drop, pick_t + post + edge_drop)
        # 5⃣ final precise window
        st.trim(pick_t - pre, pick_t + post)
        if not st:
            return None
        tr = st[0]
        tr.stats.station = sta
        tr.stats.channel = cha
        tr.stats.master_id = master_id
        return tr
    except Exception as e:
        logging.debug("Trace preprocess fail %s: %s", fpath, e)
        return None


def load_trace_for_pick(row, cfg):
    pick_t = row["pick_time"]
    ymd = pick_t.strftime("%Y%m%d")
    subdir = os.path.join(WAVEFORM_ROOT, pick_t.strftime("%Y/%m/%d"))
    pattern = f"{ymd}.*.{row['sta']}..{row['chan']}.mseed"
    matches = glob.glob(os.path.join(subdir, pattern))
    if not matches:
        return None
    return preprocess_trace(matches[0], pick_t, cfg.pre, cfg.post,
                            cfg.buffer, cfg.edge_drop,
                            cfg.band_lo, cfg.band_hi,
                            row["sta"], row["chan"], row["master_id"])

###############################################################################
# Cross‑correlation core
###############################################################################

def ncc_max(trA, trB, max_samps):
    """Return |ρ_max| and lag in samples (A vs B)."""
    cc = correlate(trA.data, trB.data, max_samps)
    norm = np.linalg.norm(trA.data) * np.linalg.norm(trB.data)
    if norm == 0:
        return 0.0, 0
    cc /= norm
    imax = np.argmax(np.abs(cc))
    return cc[imax], imax - max_samps


def crosscorr_worker(args):
    cfg, cid, net, sta, cha, pha, df_g, out_dir = args
    if len(df_g) < 2:
        return
    out_path = os.path.join(out_dir, f"cluster_{cid}_{net}_{sta}_{cha}_{pha}.csv")
    with open(out_path, "w") as f:
        f.write("oridA,oridB,station,channel,phase,lag_s,ncc\n")

    # sort by pick time for reproducibility
    df_g = df_g.sort_values("pick_time").reset_index(drop=True)
    N = len(df_g)
    n_chunks = int(np.ceil(N / cfg.chunk))
    written = set()

    def load_chunk(df_sub):
        st = []
        for _, row in df_sub.iterrows():
            tr = load_trace_for_pick(row, cfg)
            if tr is not None:
                st.append(tr)
        return st

    for i in range(n_chunks):
        i0, i1 = i * cfg.chunk, min((i + 1) * cfg.chunk, N)
        st_i = load_chunk(df_g.iloc[i0:i1])
        if not st_i:
            continue
        sr = st_i[0].stats.sampling_rate
        max_samps = int(round(cfg.max_shift * sr))

        # self‑corr inside chunk i
        for idxA in range(len(st_i) - 1):
            for idxB in range(idxA + 1, len(st_i)):
                r, lag = ncc_max(st_i[idxA], st_i[idxB], max_samps)
                if abs(r) < cfg.cc_thr:
                    continue
                a, b = st_i[idxA].stats.master_id, st_i[idxB].stats.master_id
                if (a, b) in written or (b, a) in written:
                    continue
                written.add((a, b))
                with open(out_path, "a") as f:
                    f.write(f"{a},{b},{sta},{cha},{pha},{lag/sr:.3f},{r:.3f}\n")

        # cross‑chunk i vs j>i
        for j in range(i + 1, n_chunks):
            j0, j1 = j * cfg.chunk, min((j + 1) * cfg.chunk, N)
            st_j = load_chunk(df_g.iloc[j0:j1])
            if not st_j or not st_i:
                continue
            for trA in st_i:
                for trB in st_j:
                    r, lag = ncc_max(trA, trB, max_samps)
                    if abs(r) < cfg.cc_thr:
                        continue
                    a, b = trA.stats.master_id, trB.stats.master_id
                    if (a, b) in written or (b, a) in written:
                        continue
                    written.add((a, b))
                    with open(out_path, "a") as f:
                        f.write(f"{a},{b},{sta},{cha},{pha},{lag/sr:.3f},{r:.3f}\n")
            del st_j
            gc.collect()
        del st_i
        gc.collect()
    logging.info("cluster %s %s.%s.%s done (%d pairs)", cid, sta, cha, pha, len(written))

###############################################################################
# Main routine
###############################################################################

def build_tasks(df_arr, df_ori, cfg):
    tasks = []
    for cid in sorted(df_arr["cluster_id"].unique()):
        if cid < cfg.start_cluster:
            continue
        if df_ori[df_ori["cluster_id"] == cid].shape[0] < 2:
            continue
        cdir = os.path.join(XCORR_OUTDIR, f"cluster_{cid}")
        os.makedirs(cdir, exist_ok=True)

        for (net, st, ch, ph), grp in df_arr[df_arr["cluster_id"] == cid].groupby(
            ["network", "sta", "chan", "phase"]
        ):
            if len(grp) < 2:
                continue
            tasks.append((cfg, cid, net, st, ch, ph, grp, cdir))
    return tasks


def main():
    parser = argparse.ArgumentParser(description="Cross‑correlation for GrowClust (robust version)")
    parser.add_argument("--bandpass-lo", type=float, default=BANDPASS_LO_DEF)
    parser.add_argument("--bandpass-hi", type=float, default=BANDPASS_HI_DEF)
    parser.add_argument("--pre", type=float, default=PRE_TIME_DEF)
    parser.add_argument("--post", type=float, default=POST_TIME_DEF)
    parser.add_argument("--buffer", type=float, default=BUFFER_DEF)
    parser.add_argument("--edge-drop", type=float, default=EDGE_DROP_DEF)
    parser.add_argument("--max-shift", type=float, default=MAX_SHIFT_SEC_DEF)
    parser.add_argument("--cc-threshold", type=float, default=CC_THRESHOLD_DEF)
    parser.add_argument("--chunk-size", type=int, default=CHUNK_SIZE_DEF)
    parser.add_argument("--start-cluster", type=int, default=START_CLUSTER_DEF)
    parser.add_argument("--km-radius", type=float, default=KM_RADIUS_DEF)
    parser.add_argument("--time-window", type=float, default=TIME_WINDOW_DAYS_DEF)
    parser.add_argument("--min-samples", type=int, default=MIN_SAMPLES_DEF)
    parser.add_argument("--processes", type=int, default=None)
    parser.add_argument("--log", default=LOGLEVEL_DEF, help="LOG level (DEBUG,INFO,WARNING,ERROR)")
    cfg = parser.parse_args()

    setup_logging(cfg.log)
    os.makedirs(XCORR_OUTDIR, exist_ok=True)
    logging.info("Loading & clustering …")
    df_o, df_a = load_and_cluster_dbscan(cfg.km_radius, cfg.time_window, cfg.min_samples)
    df_a = df_a.merge(df_o[["master_id", "cluster_id"]], on="master_id", how="left")
    df_a = df_a[df_a["cluster_id"] != -1].copy()
    df_a["pick_time"] = df_a["arrival_time"].apply(UTCDateTime)

    # add network if station file present
    if os.path.exists(STATIONS_PATH):
        stn = pd.read_csv(STATIONS_PATH)[["unique_id", "station", "network"]].drop_duplicates()
        stn.rename(columns={"station": "sta"}, inplace=True)
        df_a = df_a.merge(stn, on="sta", how="left")
        df_a["network"] = df_a["network"].fillna("NA")
    else:
        df_a["network"] = "NA"

    tasks = build_tasks(df_a, df_o, cfg)
    logging.info("Prepared %d tasks, spawning workers …", len(tasks))

    procs = cfg.processes or mp.cpu_count()
    with Pool(processes=procs) as pool:
        for _ in pool.imap_unordered(crosscorr_worker, tasks, chunksize=1):
            pass
    logging.info("All done. CSV files are in %s", XCORR_OUTDIR)

###############################################################################
if __name__ == "__main__":
    main()
