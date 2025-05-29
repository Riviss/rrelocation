# rrelocation

This repository contains scripts for cross‑correlating seismic waveforms and performing relative relocation using the GrowClust3D software. The code was adapted for internal use and expects data in specific directory structures.

## Repository Layout

```
src/        Python utilities and cross‑correlation helpers
scripts/    Shell and Julia run scripts
data/       Sample input/output folders
docs/       Additional documentation (see `docs/FILE_FORMATS.md`)
```

Key Python scripts now live in `src/`:

- `create_cc_and_origins.py` – cluster catalog origins with DBSCAN and generate cross‑correlation tasks.
- `create_growclust_files_after_cc.py` – build `evlist.txt`, `stlist.txt` and `xcordata.txt` files for GrowClust.
- `create_relocation_files.py` – chunked cross‑correlation workflow that limits memory use.
- `mysql_pipeline.py` – detect cluster changes, rerun relocation and upload results to a MySQL database.

The `scripts/` folder contains the Julia relocation script (`run_growclust3D.jl` and its multi‑process variant) and a helper shell script `run_relocation.sh` which loops over cluster directories.


Configuration of file paths is centralised in `src/config.py`.

## Usage

1. Generate cross‑correlation files:
   ```bash
   python src/create_cc_and_origins.py
   ```
   or use `create_relocation_files.py` for the chunked approach.

2. Create GrowClust input files for each cluster:
   ```bash
   python src/create_growclust_files_after_cc.py
   ```

3. Run the relocation for all clusters:
   ```bash
   bash scripts/run_relocation.sh
   ```

4. Upload relocated events to MySQL:
   ```bash
   python src/mysql_pipeline.py
   ```

Edit the paths in `src/config.py` to match your data layout before running the scripts.

## Installation

Create a virtual environment and install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```



## License

This project is licensed under the terms of the MIT License.  See [LICENSE](LICENSE) for details.