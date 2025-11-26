# Step 3: CTI Calculations

Purpose: extract sea surface temperature (SST) for GBIF occurrences, thin spatially dense records, compute species thermal indices (STI), and derive community temperature indices (CTI).

Primary script
- `CTI_calculation_BCODMO.R`: annotated script performing SST extraction (WorldClim WC_bio6 via sdmpredictors/terra), 50 km spatial thinning, STI summaries, and CTI computation from community matrices.

Inputs
- `../2_gbif_workflow/gbif_downloads/clean_csvs/*.csv`: cleaned GBIF occurrences.
- `../1_raw_data/outputs/pivot_all.csv`: community matrix with Year column.
- `../1_raw_data/outputs/presence_pivot_merged_sp.csv`: used for taxon alignment.

Outputs (`outputs/`)
- `species_latitudes.csv`: min/mean/max latitudes by taxon after thinning.
- `STI_results_by_taxon.csv`: STI mean/quantiles and northern latitude metrics.
- `pivot_clean.csv`: community matrix augmented with CTI metrics (mean_sti, mean_sti_2.5, mean_sti_range).

Software
- R >= 4.3 with tidyverse, terra, sf, sp, sdmpredictors, rgbif, ggrepel.

Run order
1. Confirm GBIF cleaned CSVs (Step 2) and `pivot_all.csv` (Step 1) exist.
2. Run `CTI_calculation_BCODMO.R`; check outputs in `outputs/`.
