# Usage Guide

This document provides a quick overview on how to use the demonstration script.

1. Ensure you have Python 3.8+ installed.
2. Create and activate a virtual environment (optional but recommended).
3. Run `python src/create_cc_and_origins.py --help` to view the available
   command line options. A typical invocation is:

   ```bash
   python src/create_cc_and_origins.py --start-cluster 0 --cc-threshold 0.6
   ```

   Adjust `--start-cluster` or `--cc-threshold` as needed and use `--processes`
   to control the number of worker processes.

As development continues, this guide will be expanded with more detailed
instructions and explanations.