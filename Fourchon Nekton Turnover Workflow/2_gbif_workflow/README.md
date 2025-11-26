# Step 2: GBIF Downloads and Cleaning

## Abstract
This step assembles GBIF occurrence records for all taxa detected in the Port Fourchon nekton surveys, cleans spatial artifacts, and delivers per-taxon CSVs ready for thermal niche extraction. The workflow expands harmonized taxon names into GBIF queries, submits authenticated bulk downloads, filters out problematic coordinates and land points, and documents all DOI citations needed for downstream analyses and reporting.

Purpose: fetch occurrence data for target taxa, apply coordinate/data-quality filters, and save cleaned per-taxon CSVs with citation logs.

Primary script
- `gbif_download_BCODMO.R`: annotated script that expands pooled taxa, submits GBIF downloads, cleans coordinates (CoordinateCleaner), and logs citations/DOIs.

Inputs
- `../1_raw_data/outputs/presence_pivot_merged_sp.csv`: taxon list used to build download targets.
- Environment variables: `GBIF_USER`, `GBIF_EMAIL`, `GBIF_PWD` for GBIF API access.

Outputs
- `gbif_downloads/clean_csvs/<taxon>_clean.csv`: cleaned occurrences per taxon.
- `gbif_citations.txt`: citations/DOIs for all downloads.


## Software 
    R version 4.5.1 (2025-06-13 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26100)

    Matrix products: default
    LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8
    [2] LC_CTYPE=English_United States.utf8
    [3] LC_MONETARY=English_United States.utf8
    [4] LC_NUMERIC=C
    [5] LC_TIME=English_United States.utf8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] CoordinateCleaner_3.0.1 sdmpredictors_0.2.15    sp_2.2-0
    [4] terra_1.8-60            lubridate_1.9.4         forcats_1.0.1
    [7] stringr_1.5.2           dplyr_1.1.4             purrr_1.1.0
    [10] readr_2.1.5             tidyr_1.3.1             tibble_3.3.0
    [13] ggplot2_4.0.0           tidyverse_2.0.0         rgbif_3.8.3

    loaded via a namespace (and not attached):
    [1] gtable_0.3.6        jsonlite_2.0.0      compiler_4.5.1
    [4] tidyselect_1.2.1    Rcpp_1.1.0          geosphere_1.5-20
    [7] xml2_1.4.0          scales_1.4.0        lattice_0.22-7
    [10] R6_2.6.1            plyr_1.8.9          generics_0.1.4
    [13] oai_0.4.0           tzdb_0.5.0          pillar_1.11.1
    [16] RColorBrewer_1.1-3  rlang_1.1.6         stringi_1.8.7
    [19] S7_0.2.0            rnaturalearth_1.1.0 lazyeval_0.2.2
    [22] timechange_0.3.0    cli_3.6.5           withr_3.0.2
    [25] magrittr_2.0.3      grid_4.5.1          hms_1.1.4
    [28] lifecycle_1.0.4     vctrs_0.6.5         glue_1.8.0
    [31] data.table_1.17.8   raster_3.6-32       farver_2.1.2
    [34] whisker_0.4.1       codetools_0.2-20    httr_1.4.7
    [37] tools_4.5.1         pkgconfig_2.0.3

Run order
1. Ensure Step 1 outputs are present.
2. Set GBIF credentials in the environment.
3. Run `gbif_download_BCODMO.R`; verify cleaned CSVs and `gbif_citations.txt`.
