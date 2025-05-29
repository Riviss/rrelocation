# rrelocation

This repository contains scripts for cross‑correlating seismic waveforms and performing relative relocation using the GrowClust3D software. The code was adapted for internal use and expects data in specific directory structures.

## Overview of Scripts

- **create_cc_and_origins.py** – Generates cross‑correlation tasks by clustering catalog origins using DBSCAN and then correlating waveform pairs. Outputs files under `xcorr_output/`.
- **create_growclust_filles_after_cc.py** – Converts cross‑correlation results into the `evlist.txt`, `stlist.txt`, and `xcordata.txt` files required by GrowClust.
- **create_relocation_files.py** – Alternative cross‑correlation workflow that processes waveforms in chunks to limit memory use.
- **run_growclust3D.jl** and **run_growclust3D-MP.jl** – Julia scripts that run the GrowClust3D relocation algorithm in single- and multi‑process modes.
- **run_relocation.sh** – Convenience shell script that loops over cluster directories and calls the Julia relocation script.

## Usage

1. Generate cross‑correlation files:
   ```bash
   python create_cc_and_origins.py
   ```
   or use `create_relocation_files.py` for the chunked approach.

2. Create GrowClust input files for each cluster:
   ```bash
   python create_growclust_filles_after_cc.py
   ```

3. Run the relocation for all clusters:
   ```bash
   bash run_relocation.sh
   ```

Edit the paths at the top of each script to match your data layout before running.

## License

This project is licensed under the terms of the MIT License.  See [LICENSE](LICENSE) for details.