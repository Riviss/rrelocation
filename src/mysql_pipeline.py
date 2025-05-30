#!/usr/bin/env python3
"""MySQL relocation pipeline.

This utility clusters events, checks for changes and runs relocation
for affected clusters. Relocated events are then uploaded to a MySQL
server in the ``earthquakes`` database.

The script expects the usual input CSV files configured in ``config.py``
(the ``ORIGINS_PATH`` and ``ARRIVALS_PATH``). Cluster information is
stored in a table named ``event_clusters`` and relocation results are
written to ``relocated_events``.

Edit the database credentials as needed before running. Connections are
established using SQLAlchemy and default to ``root@localhost`` with an empty
password.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from datetime import datetime
from typing import Dict, Iterable, List, Set
import multiprocessing
from multiprocessing import Pool

import pandas as pd

from sqlalchemy import create_engine, text
from obspy import UTCDateTime

from config import (
    AGG_CAT,
    AGG_CLUST,
    AGG_LOG,
    GROWCLUST_INP,
    JULIA_SCRIPT,
    XCORR_OUTDIR,
    USE_DATABASE,
    DB_HOST,
    DB_USER,
    DB_PASSWORD,
    DB_NAME,
    STATIONS_CSV,
    CLUSTER_ROOT,
    get_db_engine,
)

from create_cc_and_origins import load_and_cluster_dbscan, crosscorr_group
from create_growclust_files_after_cc import (
    create_masterid_map,
    create_station_mapping,
    write_evlist,
    write_stlist,
    write_xcordata,
    load_station_info,
)
###############################################################################
# Database helpers
###############################################################################

def get_engine() -> "Engine":
    """Return a SQLAlchemy engine using config credentials."""

    if USE_DATABASE:
        return get_db_engine()
    # Fallback for backwards compatibility
    if DB_PASSWORD:
        auth = f"{DB_USER}:{DB_PASSWORD}"
    else:
        auth = DB_USER
        url = f"mysql+mysqlconnector://{auth}@{DB_HOST}/{DB_NAME}"
        return create_engine(url)


def ensure_tables(engine) -> None:
    """Create required tables if they do not exist."""
    with engine.begin() as conn:
        conn.execute(
            text(
                """
                CREATE TABLE IF NOT EXISTS event_clusters (
                    master_id    VARCHAR(64) PRIMARY KEY,
                    cluster_id   INT,
                    lat          DOUBLE,
                    lon          DOUBLE,
                    depth        DOUBLE,
                    mag          DOUBLE,
                    last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                                     ON UPDATE CURRENT_TIMESTAMP
                )
                """
            )
        )
        conn.execute(
            text(
                """
                CREATE TABLE IF NOT EXISTS relocated_events (
                    master_id  INT PRIMARY KEY,
                    datetime   DATETIME,
                    lat        DOUBLE,
                    lon        DOUBLE,
                    depth      DOUBLE,
                    mag        DOUBLE,
                    rmsP       DOUBLE,
                    rmsS       DOUBLE,
                    eh         DOUBLE,
                    ez         DOUBLE,
                    et         DOUBLE
                )
                """
            )
        )


def fetch_previous_clusters(engine) -> Dict[str, Dict[str, float]]:
    with engine.connect() as conn:
        rows = conn.execute(text("SELECT * FROM event_clusters")).mappings().all()
    prev: Dict[str, Dict[str, float]] = {}
    for r in rows:
        prev[str(r["master_id"])] = {
            "cluster_id": int(r["cluster_id"]),
            "lat": float(r["lat"]),
            "lon": float(r["lon"]),
            "depth": float(r["depth"]),
        }
    return prev


def update_event_clusters(engine, df: pd.DataFrame) -> None:
    q = text(
        "REPLACE INTO event_clusters (master_id, cluster_id, lat, lon, depth, mag) "
        "VALUES (:master_id, :cluster_id, :lat, :lon, :depth, :mag)"
    )
    with engine.begin() as conn:
        for _, row in df.iterrows():
            conn.execute(
                q,
                {
                    "master_id": str(row["master_id"]),
                    "cluster_id": int(row["cluster_id"]),
                    "lat": float(row["lat"]),
                    "lon": float(row["lon"]),
                    "depth": float(row["depth"]),
                    "mag": float(row.get("mag", float("nan"))),
                },
            )


def upload_relocated_events(engine, cat_file: str) -> None:
    if not os.path.exists(cat_file):
        print(f"No relocation catalogue found at {cat_file}")
        return
    cols = [
        "year",
        "month",
        "day",
        "hour",
        "minute",
        "second",
        "master_id",
        "lat",
        "lon",
        "depth",
        "mag",
        "qID",
        "cID",
        "nbranch",
        "qnpair",
        "qndiffP",
        "qndiffs",
        "rmsP",
        "rmsS",
        "eh",
        "ez",
        "et",
        "latC",
        "lonC",
        "depC",
    ]
    df = pd.read_csv(cat_file, delim_whitespace=True, header=None, names=cols)
    q = text(
        "REPLACE INTO relocated_events "
        "(master_id, datetime, lat, lon, depth, mag, rmsP, rmsS, eh, ez, et) "
        "VALUES (:master_id, :datetime, :lat, :lon, :depth, :mag, :rmsP, :rmsS, :eh, :ez, :et)"
    )
    with engine.begin() as conn:
        for _, r in df.iterrows():
            dt = datetime(
                int(r["year"]),
                int(r["month"]),
                int(r["day"]),
                int(r["hour"]),
                int(r["minute"]),
                int(float(r["second"])),
            )
            conn.execute(
                q,
                {
                    "master_id": int(r["master_id"]),
                    "datetime": dt,
                    "lat": float(r["lat"]),
                    "lon": float(r["lon"]),
                    "depth": float(r["depth"]),
                    "mag": float(r["mag"]),
                    "rmsP": float(r["rmsP"]),
                    "rmsS": float(r["rmsS"]),
                    "eh": float(r["eh"]),
                    "ez": float(r["ez"]),
                    "et": float(r["et"]),
                },
            )
    print(f"Uploaded {len(df)} relocated events")


###############################################################################
# Cross-correlation and file builders
###############################################################################

def run_crosscorr_for_clusters(
    df_origin: pd.DataFrame,
    df_arrival: pd.DataFrame,
    cluster_ids: Iterable[int],
    processes: int | None = None,
) -> None:
    """Run cross-correlation for selected clusters only."""

    arr = pd.merge(
        df_arrival,
        df_origin[["master_id", "cluster_id"]],
        on="master_id",
        how="left",
    )
    arr = arr[arr["cluster_id"] != -1].copy()
    arr["pick_time"] = arr["arrival_time"].apply(lambda x: UTCDateTime(x))

    if os.path.exists(STATIONS_CSV):
        df_st = pd.read_csv(STATIONS_CSV)[["unique_id", "station", "network"]].drop_duplicates()
        df_st.rename(columns={"station": "sta"}, inplace=True)
        arr = pd.merge(arr, df_st, on="sta", how="left")
        arr["network"] = arr["network"].fillna("NA")
    else:
        arr["network"] = "NA"

    cpu_count = processes or multiprocessing.cpu_count()

    tasks = []
    for cid in cluster_ids:
        dfc = arr[arr["cluster_id"] == cid]
        if dfc.empty:
            continue
        n_events = df_origin[df_origin["cluster_id"] == cid].shape[0]
        if n_events < 2:
            continue

        cdir = os.path.join(XCORR_OUTDIR, f"cluster_{cid}")
        os.makedirs(cdir, exist_ok=True)

        for (net, st, ch, ph), df_g in dfc.groupby(["network", "sta", "chan", "phase"]):
            if len(df_g) < 2:
                continue
            tasks.append(
                {
                    "cluster_id": cid,
                    "network": net,
                    "sta": st,
                    "chan": ch,
                    "phase": ph,
                    "df_group": df_g,
                    "cluster_dir": cdir,
                }
            )

    if not tasks:
        print("No cross-correlation tasks for changed clusters")
        return

    with Pool(processes=cpu_count) as pool:
        results = pool.map(crosscorr_group, tasks)

    for msg in results:
        print("Worker =>", msg)


def build_growclust_files_for_clusters(
    df_origin: pd.DataFrame, cluster_ids: Iterable[int]
) -> None:
    """Create GrowClust input files for selected clusters."""

    df_events = df_origin[["master_id", "lat", "lon", "depth", "datetime", "cluster_id"]].copy()
    df_stations = load_station_info(STATIONS_CSV)
    master_map = create_masterid_map(df_events)
    station_map = create_station_mapping(df_stations)

    for cid in cluster_ids:
        cdir = os.path.join(CLUSTER_ROOT, f"cluster_{cid}")
        os.makedirs(cdir, exist_ok=True)
        dfc = df_events[df_events["cluster_id"] == cid].copy()
        if dfc.empty:
            continue
        write_evlist(dfc, cid, cdir, master_map)
        write_stlist(df_stations, cid, cdir)
        write_xcordata(cid, cdir, master_map, station_map)




###############################################################################
# Relocation runner
###############################################################################

def run_relocation_for_clusters(cluster_ids: Iterable[int]) -> None:
    os.makedirs("data/in", exist_ok=True)
    os.makedirs("data/out", exist_ok=True)
    os.makedirs(os.path.dirname(AGG_CAT), exist_ok=True)

    # Empty aggregated outputs
    open(AGG_CAT, "w").close()
    open(AGG_CLUST, "w").close()
    open(AGG_LOG, "w").close()

    for cid in cluster_ids:
        cdir = os.path.join(XCORR_OUTDIR, f"cluster_{cid}")
        if not os.path.isdir(cdir):
            print(f"Cluster dir {cdir} missing; skipping")
            continue

        shutil.copy(os.path.join(cdir, "evlist.txt"), "data/in/evlist.txt")
        shutil.copy(os.path.join(cdir, "stlist.txt"), "data/in/stlist.txt")
        shutil.copy(os.path.join(cdir, "xcordata.txt"), "data/in/xcordata.txt")

        for f in [
            "data/out/out.nllgrid3D.cat",
            "data/out/out.nllgrid3D.clust",
            "data/out/out.nllgrid3D.log",
        ]:
            if os.path.exists(f):
                os.remove(f)

        subprocess.run(["julia", JULIA_SCRIPT, GROWCLUST_INP], check=True)

        with open(AGG_CAT, "a") as acat, open("data/out/out.nllgrid3D.cat") as src:
            acat.write(src.read())
        with open(AGG_CLUST, "a") as acl, open("data/out/out.nllgrid3D.clust") as src:
            acl.write(src.read())
        with open(AGG_LOG, "a") as alog, open("data/out/out.nllgrid3D.log") as src:
            alog.write(src.read())


###############################################################################
# Cluster change detection
###############################################################################

def detect_changed_clusters(
    df_origin: pd.DataFrame, prev: Dict[str, Dict[str, float]]
) -> List[int]:
    changed: Set[int] = set()
    seen_prev: Set[str] = set()
    for _, row in df_origin.iterrows():
        eid = str(row["master_id"])
        cid = int(row["cluster_id"])
        lat = float(row["lat"])
        lon = float(row["lon"])
        dep = float(row["depth"])
        info = prev.get(eid)
        if info is None:
            changed.add(cid)
            continue
        seen_prev.add(eid)
        loc_changed = (
            abs(info["lat"] - lat) > 1e-4
            or abs(info["lon"] - lon) > 1e-4
            or abs(info["depth"] - dep) > 0.1
        )
        if info["cluster_id"] != cid or loc_changed:
            changed.add(cid)

    removed_ids = set(prev.keys()) - seen_prev
    for rid in removed_ids:
        changed.add(int(prev[rid]["cluster_id"]))

    return sorted(changed)


###############################################################################
# Main entry
###############################################################################

def main() -> None:
    engine = get_engine()
    ensure_tables(engine)

    origin_df, arrival_df = load_and_cluster_dbscan()
    prev = fetch_previous_clusters(engine)
    changed_clusters = detect_changed_clusters(origin_df, prev)

    if not changed_clusters:
        print("No cluster changes detected")
        return

    print(f"Processing clusters: {changed_clusters}")

    run_crosscorr_for_clusters(origin_df, arrival_df, changed_clusters)
    build_growclust_files_for_clusters(origin_df, changed_clusters)
    run_relocation_for_clusters(changed_clusters)
    upload_relocated_events(engine, AGG_CAT)
    update_event_clusters(engine, origin_df)


if __name__ == "__main__":
    main()