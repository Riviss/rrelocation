"""Centralized file path configuration for the relocation workflow."""
import os

# Input data

# Use database instead of CSV files
USE_DATABASE = False

# Database connection
DB_HOST = "localhost"
DB_USER = "root"
DB_PASSWORD = None
DB_NAME = "earthquakes"
ORIGINS_TABLE = "master_origin_3D"
ARRIVALS_TABLE = "master_arrival_3D"


def get_db_engine():
    """Return SQLAlchemy engine using the configured credentials."""
    if DB_PASSWORD:
        auth = f"{DB_USER}:{DB_PASSWORD}"
    else:
        auth = DB_USER
    url = f"mysql+mysqlconnector://{auth}@{DB_HOST}/{DB_NAME}"
    return create_engine(url)

ORIGINS_PATH = "/home/pgcseiscomp/Documents/projects/velocitymodel_to_traveltimegrid/events/KSMMA/origins.csv"
ARRIVALS_PATH = "/home/pgcseiscomp/Documents/projects/velocitymodel_to_traveltimegrid/events/KSMMA/arrivals.csv"
STATIONS_PATH = "/home/pgcseiscomp/Documents/seismic_process/velocity_model/nll/stations/stations.csv"
WAVEFORM_ROOT = "/home/pgcseiscomp/antelope/wfs"

# Cross-correlation output directory
XCORR_OUTDIR = "./xcorr_output"
os.makedirs(XCORR_OUTDIR, exist_ok=True)

# Additional files used after cross-correlation
CLUSTER_ROOT = XCORR_OUTDIR
EVENTS_CSV = os.path.expanduser("~/Documents/seismic_process/relocation_3D/nebc/xcorr_output/clustered_origins.csv")
STATIONS_CSV = STATIONS_PATH  # used by create_growclust_filles_after_cc

# GrowClust related
JULIA_SCRIPT = "/home/pgcseiscomp/Documents/seismic_process/relocation_3D/nebc/run_growclust3D.jl"
GROWCLUST_INP = "dawson.nllgrid3D.inp"

# Aggregated output locations used by run_relocation.sh
AGG_CAT = "data/aggregated/out.nllgrid3D.cat"
AGG_CLUST = "data/aggregated/out.nllgrid3D.clust"
AGG_LOG = "data/aggregated/out.nllgrid3D.log"