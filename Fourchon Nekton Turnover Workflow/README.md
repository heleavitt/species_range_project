# Data for  *Resilient Nekton Composition in the Face of Climate-Driven Foundation Species Shifts* 

This package contains data and scripts use to run analyses and produce figures for publication Leavitt, H; Thomas, A; Doerr, J; Johnson, D; Nelson, J. (In press) Resilient Nekton Composition in the Face of Climate-Driven Foundation Species Shifts. Ecology. Accepted 2025-11-14

## Software Versions: 
### R version 4.5.1 (2025-06-13 ucrt)
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
    [1] sdmpredictors_0.2.15

    loaded via a namespace (and not attached):
    [1] compiler_4.5.1   cli_3.6.5        tools_4.5.1      sp_2.2-0        
    [5] Rcpp_1.1.0       raster_3.6-32    codetools_0.2-20 grid_4.5.1
    [9] jsonlite_2.0.0   rlang_1.1.6      lattice_0.22-7   terra_1.8-60

###  Python version: 3.13.0 | packaged by Anaconda, Inc. | (main, Oct  7 2024, 21:21:52) [MSC v.1929 64 bit (AMD64)]
    cartopy==0.25.0
    matplotlib==3.10.7
    numpy==2.3.5
    pandas==2.3.3
    pandasgui==0.2.15
    Requests==2.32.5
    seaborn==0.13.2


## The workflow proceeds in four numbered steps:
- `1_raw_data/`: Merge raw drop sample datasets (2006, 2016, 2022, 2023) and produce community/presence tables.
- `2_gbif_workflow/`: Download GBIF occurrences for target taxa, clean coordinates, and save per-taxon CSVs plus citations.
- `3_CTI_calculations/`: Compute species thermal indices (STI) and community temperature indices (CTI) using GBIF and Wordclim climate data extractions.
- `4_species_of_interest/`: Apply sequential filters to identify co-migrator candidates and plot CTI against winter temperature.

### Run order
1. In `1_raw_data/`, run `dataset_merge_and_process_BCODMO.ipynb`

   **Inputs**: 
    -   ` inputs/raw_2006_fish.csv`, inputs/raw_2006_invert.csv, inputs/raw_2006_peneid.csv, inputs/2006_site_data.csv - Data from 2006 was split up by clades with separate information on each site. Data is in long format
    -   ` inputs/raw_2016.csv` - raw wide-format data from 2016 drop sampling study by Doerr et al. 
    -   ` inputs/raw_2022.csv`, inputs/raw_2023.csv - raw wide-format data from 2022 and 2023 drop sampling studies 
    -   ` inputs/final_aphia_codex_edited.csv` - hand-generated list of species observed in each study harmonized to current valid taxonomic names as of 2025 along with aphia IDs used in WoRMS database. 
    
    
    **Functions**: 
    script wrangles different datasets into a standardized format for analysis. Also compiles some species into genera or families based on the lowest common resolution of all datasets or changes to taxonomix classifications over time. 
    
    
    **Outputs**: 
    -    `pivot_all.csv` - pivot table of includeing data from all three studies wioth sites as rows as rows and species as columns.     
    -    `presence_summary.csv` and `presence_summary2.csv` - table with taxomic groups in each study as rows (ex. Minuca spp. in 2006 is one row, Minuca spp. in 2022_2023 is another). Includes column of propotion of sites each species is observed in that study. Presence_summary2 includes harmonized taxonmic groups, while presency_summary includes origional taxonomic groups. 

    -    `presence_pivot_merged.csv` - table with species as rows, and studies as columns. Values are proportion of sites each species was observed at in that study. 

2. In `2_gbif_workflow/`, run `gbif_download_BCODMO.R` 

    **Inputs**: `presence_pivot_merged.csv`


    **Functions**: Query GBIF for observations of all species oberved throughout all three studies. Clean observations for use in CTI calculations later on. 


    **Outputs**: gbif_download folder includes raw GBIF downloads for each species as .zip files. clean_csvs folder includes .csv files for each species post cleanup. Text file of doi citations for each species is also included in this folder. 

3. In `3_CTI_calculations/`, run `CTI_calculation_BCODMO.R` 
    
    **Inputs**: 
    -    `2_gbif_workflow/clean_csvs/*_clean.csv` - cleaned GBIF occurrences per taxon.
    -    `pivot_all.csv`- community matrix with SampleID + Year.
    -    `species_presence_2006_2016_2022.csv`- presence metadata used in joins.
    
    
    **Function**:  Read cleaned GBIF occurrence records for focal taxa.
    - Extract mean sea surface temperature (SST) for each occurrence.
    - Thin spatially dense records to reduce sampling bias.
    - Derive per-species thermal metrics (STI and confidence bounds).
    - Combine with community matrices to compute CTI-style site metrics.
    
    **Outputs**:     
    - outputs/species_latitudes.csv         : latitudinal summaries for each taxon.
    - outputs/STI_results_by_taxon.csv      : STI statistics per taxon after thinning.
    - outputs/pivot_clean.csv               : community matrix with CTI metrics appended.

4. In `4_species_of_interest/`, run `ID_species_of_interest_BCODMO.R` 

    **inputs**
    - `pivot_clean.csv`              : community matrix with CTI metrics attached.
    - `STI_results_by_taxon.csv`     : species thermal indices from GBIF extractions.
    - `species_presence_2006_2016_2022.csv`   : coarse presence by year.
    - `avg_winter_temp_yearly.csv`                     : winter temperature time series (Port Fourchon).

    **Function**: Run species pool through a series of filters to identify species that fit description of co-migrator with mangroves.

    **Outputs**:
    - outputs/filter1_species.csv         : species that are increasing in prevalence over time or have only been cited in later studies
    - outputs/filter2_species.csv         : species that have a sufficient number of GBIF observations in the data to be confident of CTI and other metrics
    - outputs/filter3_species.csv         : species that show 97.5% of observations in climates warmer than port fourchon (warm adapted species at the edge of their range)
    - outputs/filter4_species.csv         : near northern edge of observed range (99th percentile lat < 30.7)
    - outputs/filter_final.csv            : species passing all filters and have a low chance of false absences
    - outputs/plots/species_temp_plot.png : CTI vs winter temperature figure.


