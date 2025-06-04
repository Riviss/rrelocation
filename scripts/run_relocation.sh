#!/bin/bash

# Create necessary directories if they don't exist.
mkdir -p data/in
mkdir -p data/out
mkdir -p data/aggregated

# Aggregated output files (final collected outputs)
agg_cat="data/aggregated/out.nllgrid3D.cat"
agg_clust="data/aggregated/out.nllgrid3D.clust"
agg_log="data/aggregated/out.nllgrid3D.log"

# Initialize the aggregated files.
> "$agg_cat"
> "$agg_clust"
> "$agg_log"

# Determine repository root based on this script location
script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"

# Loop through each cluster folder in xcorr_output.
for cluster_folder in xcorr_output/cluster_*; do
  [ -d "$cluster_folder" ] || continue

  # Clear previous input files.
  rm -f data/in/*

  # Copy the required files from the current cluster folder to data/in.
  cp "$cluster_folder/evlist.txt" data/in/evlist.txt
  cp "$cluster_folder/stlist.txt" data/in/stlist.txt
  cp "$cluster_folder/xcordata.txt" data/in/xcordata.txt

  # Remove previous Julia outputs if present.
  rm -f data/out/out.nllgrid3D.cat data/out/out.nllgrid3D.clust data/out/out.nllgrid3D.log

  # Run the Julia script using paths relative to the repository
  julia "$repo_root/scripts/run_growclust3D.jl" \
        "$repo_root/dawson.nllgrid3D.inp"

  # Use the folder name (e.g., cluster_0) as a marker.
  cluster_name=$(basename "$cluster_folder")
  #echo "### $cluster_name" >> "$agg_cat"
  echo "### $cluster_name" >> "$agg_clust"
  echo "### $cluster_name" >> "$agg_log"

  # Append the current run's Julia output into aggregated files.
  cat data/out/out.nllgrid3D.cat >> "$agg_cat"
  cat data/out/out.nllgrid3D.clust >> "$agg_clust"
  cat data/out/out.nllgrid3D.log >> "$agg_log"
done

