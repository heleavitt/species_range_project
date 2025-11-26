# Step 1: Raw Data Merge and Presence Tables

Purpose: consolidate seine catch datasets (2006, 2016, 2022, 2023) and produce community matrices and presence summaries for downstream CTI analyses.

Primary script/notebook
- `dataset_merge_and_process_BCODMO.ipynb`: annotated Jupyter notebook that loads raw CSVs, standardizes taxonomic names, and outputs presence/community products.

Inputs (`inputs/`)
- `raw_2006_fish.csv`, `raw_2006_invert.csv`, `raw_2006_peneid.csv`, `2006_site_data.csv` (+ parameter descriptions).
- `raw_2016.csv` (+ parameter description).
- `raw_2022.csv`, `raw_2023.csv` (+ parameter descriptions).
- `final_aphia_codex_edited.csv` (taxonomic lookup across survey years).
- `species_presence_2006_2016_2022.csv` (coarse presence metadata).

Outputs (`outputs/`)
- `pivot_all.csv`: community matrix (SampleID x Taxon) with Year column.
- `presence_summary.csv`, `presence_pivot.csv`: per-taxon site-level presence proportions.
- `presence_summary2.csv`, `presence_pivot_merged_sp.csv`: presence after taxonomic harmonization for CTI use.

Software
- Python >= 3.9; pandas. The notebook sets the working directory to the repository root.

Run order
1. Open and run `dataset_merge_and_process_BCODMO.ipynb`.
2. Confirm outputs in `outputs/` before proceeding to GBIF downloads (Step 2).
