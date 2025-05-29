# rrelocation

This repository contains scripts for cross‑correlating seismic waveforms and performing relative relocation using the `GrowClust3D.jl` package. The code was adapted for internal use and expects data in specific directory structures. `GrowClust3D.jl` is the Julia implementation of the original GrowClust algorithm.

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

The `scripts/` folder contains the Julia relocation script (`run_growclust3D.jl` and its multi‑process variant) and a helper shell script `run_relocation.sh` which loops over cluster directories. These scripts call the `GrowClust3D.jl` package, a Julia port of the original GrowClust relocation algorithm.

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

We recommend installing the Python dependencies inside an Anaconda environment:

```bash
conda create -n rreloc python=3.8
conda activate rreloc
pip install -r requirements.txt
```

## Running the Julia relocation script

The relocation routine itself is implemented in Julia using the
`GrowClust3D.jl` package.  You will need a working Julia installation
and the package available in your Julia environment.  A typical setup
looks like:

```bash
# install Julia packages (from the Julia REPL)
import Pkg; Pkg.add("GrowClust3D")

# run GrowClust with the provided control file
julia scripts/run_growclust3D.jl dawson.nllgrid3D.inp
```

Edit the `dawson.nllgrid3D.inp` control file or supply your own path if
needed.  When running batch relocations `scripts/run_relocation.sh`
will execute the above command for each cluster directory.


## License

This project is licensed under the terms of the MIT License.  See [LICENSE](LICENSE) for details.