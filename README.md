# rrelocation

This repository contains scripts for cross‑correlating seismic waveforms and performing relative relocation using the `GrowClust3D.jl` package. The code was adapted for internal use and expects data in specific directory structures. `GrowClust3D.jl` is the Julia implementation of the original GrowClust algorithm.

## Repository Layout

```
src/        Python utilities and cross‑correlation helpers
scripts/    Shell and Julia run scripts
data/       Sample input/output folders
docs/       Additional documentation
```

Key Python scripts now live in `src/`:

- `create_cc_and_origins.py` – cluster catalog origins with DBSCAN and generate cross‑correlation tasks.
- `create_growclust_files_after_cc.py` – build `evlist.txt`, `stlist.txt` and `xcordata.txt` files for GrowClust.
- `create_relocation_files.py` – chunked cross‑correlation workflow that limits memory use.

The `scripts/` folder contains the Julia relocation script (`run_growclust3D.jl` and its multi‑process variant) and a helper shell script `run_relocation.sh` which loops over cluster directories. These scripts call the `GrowClust3D.jl` package, a Julia port of the original GrowClust relocation algorithm.

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

Edit the paths at the top of each script to match your data layout before running.

## Installation

We recommend installing the Python dependencies inside an Anaconda environment:

```bash
conda create -n rreloc python=3.8
conda activate rreloc
pip install -r requirements.txt
```



## License

This project is licensed under the terms of the MIT License.  See [LICENSE](LICENSE) for details.