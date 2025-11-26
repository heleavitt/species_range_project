# Step 2: GBIF Downloads and Cleaning

Purpose: fetch occurrence data for target taxa, apply coordinate/data-quality filters, and save cleaned per-taxon CSVs with citation logs.

Primary script
- `gbif_download_BCODMO.R`: annotated script that expands pooled taxa, submits GBIF downloads, cleans coordinates (CoordinateCleaner), and logs citations/DOIs.

Inputs
- `../1_raw_data/outputs/presence_pivot_merged_sp.csv`: taxon list used to build download targets.
- Environment variables: `GBIF_USER`, `GBIF_EMAIL`, `GBIF_PWD` for GBIF API access.

Outputs
- `gbif_downloads/clean_csvs/<taxon>_clean.csv`: cleaned occurrences per taxon.
- `gbif_citations.txt`: citations/DOIs for all downloads.

Software
- R >= 4.3 with rgbif, tidyverse, terra, sp, sdmpredictors, CoordinateCleaner.

Run order
1. Ensure Step 1 outputs are present.
2. Set GBIF credentials in the environment.
3. Run `gbif_download_BCODMO.R`; verify cleaned CSVs and `gbif_citations.txt`.
