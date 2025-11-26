# Step 3: CTI Calculations

## Abstract
This step links cleaned GBIF occurrences to climatology to quantify species- and community-scale thermal affinities. Occurrence points are spatially thinned to reduce sampling bias, paired with WorldClim SST (bio6) to derive per-taxon species thermal indices and latitudinal envelopes, then joined to multi-year community matrices to compute CTI metrics (means and confidence bounds) for each sample.

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
    [1] ggrepel_0.9.6           sf_1.0-21               CoordinateCleaner_3.0.1
    [4] sdmpredictors_0.2.15    sp_2.2-0                terra_1.8-60
    [7] lubridate_1.9.4         forcats_1.0.1           stringr_1.5.2
    [10] dplyr_1.1.4             purrr_1.1.0             readr_2.1.5
    [13] tidyr_1.3.1             tibble_3.3.0            ggplot2_4.0.0
    [16] tidyverse_2.0.0         rgbif_3.8.3

    loaded via a namespace (and not attached):
    [1] generics_0.1.4      class_7.3-23        xml2_1.4.0
    [4] KernSmooth_2.23-26  stringi_1.8.7       lattice_0.22-7
    [7] hms_1.1.4           magrittr_2.0.3      grid_4.5.1
    [10] timechange_0.3.0    RColorBrewer_1.1-3  rnaturalearth_1.1.0
    [13] plyr_1.8.9          jsonlite_2.0.0      whisker_0.4.1
    [16] e1071_1.7-16        DBI_1.2.3           httr_1.4.7
    [19] scales_1.4.0        oai_0.4.0           codetools_0.2-20
    [22] lazyeval_0.2.2      cli_3.6.5           rlang_1.1.6
    [25] units_0.8-7         withr_3.0.2         tools_4.5.1
    [28] raster_3.6-32       tzdb_0.5.0          geosphere_1.5-20
    [31] vctrs_0.6.5         R6_2.6.1            proxy_0.4-27
    [34] classInt_0.4-11     lifecycle_1.0.4     pkgconfig_2.0.3
    [37] pillar_1.11.1       gtable_0.3.6        data.table_1.17.8
    [40] glue_1.8.0          Rcpp_1.1.0          tidyselect_1.2.1
    [43] farver_2.1.2        compiler_4.5.1      S7_0.2.0
    >
Run order
1. Confirm GBIF cleaned CSVs (Step 2) and `pivot_all.csv` (Step 1) exist.
2. Run `CTI_calculation_BCODMO.R`; check outputs in `outputs/`.
