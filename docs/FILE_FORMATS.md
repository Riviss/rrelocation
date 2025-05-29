# File Formats

This document describes the expected structure of the various input and output files referenced by the Python scripts.

## origins.csv
A comma separated file with at least the columns:

- `master_id` – unique identifier for each event
- `event_datetime` (renamed to `datetime` by the scripts)
- `latitude` (`lat`)
- `longitude` (`lon`)
- `depth_m` (`depth`)

## arrivals.csv
Contains pick information. Required columns include:

- `master_id` – event identifier matching `origins.csv`
- `sta` – station code
- `phase` – phase name
- `datetime` or `arrival_time` – arrival timestamp

## stations.csv
Station metadata used for mapping station codes. Expected columns are:

- `station` or `sta` – original station code
- `unique_id` – name used in relocation
- `latitude`/`longitude` (`lat`/`lon`)
- `elevation` (`elev`)

## Cross-correlation output
Text files named `xcorr_cluster_<id>_*.txt` produced by the cross-correlation
workflow. Each line has the form:

```
EVENT_A EVENT_B STATION TDIF CC PHASE
```

Example:

```
WORKING_cat_A_5367 WORKING_cat_A_1234 BCH1A -0.100 0.965 S
```

Header lines beginning with `#` are ignored.

## GrowClust input files
`create_growclust_filles_after_cc.py` generates three files in each cluster directory:

- **evlist.txt** – one event per line following
  `YYYY MM DD HH MM SS.SSS LAT LON DEPTH PLACEHOLDERS EVENT_ID`.
- **stlist.txt** – station list in the format `STA LAT LON ELEV`.
- **xcordata.txt** – aggregated cross-correlation data using numeric event IDs.

Waveform data are expected under
`/home/pgcseiscomp/antelope/wfs/YYYY/MM/DD/YYYYMMDD.*.STA..CHAN.mseed`.