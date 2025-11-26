# Step 4: Species of Interest and CTI Plotting

## Abstract
This step screens the Port Fourchon nekton assemblage to pinpoint potential mangrove co‑migrators. It chains abundance/presence trends across survey years with GBIF sampling sufficiency, thermal niche position, and range-edge proximity, then flags taxa with low false-absence risk. The pipeline outputs species retained at each filter stage and visualizes CTI trajectories against local winter temperatures to highlight warm-adapted taxa advancing into the study area.

Purpose: identify co-migrator candidates via sequential filters (abundance trends, GBIF sufficiency, thermal limits, range position, false-absence probability) and plot CTI alongside winter temperatures.

Primary script
- `ID_species_of_interest_BCODMO.R`: annotated script applying filters and generating CTI vs winter temperature plot.

Inputs
- `../3_CTI_calculations/outputs/pivot_clean.csv`: community matrix with CTI metrics.
- `../3_CTI_calculations/outputs/STI_results_by_taxon.csv`: species thermal indices.
- `../1_raw_data/inputs/species_presence_2006_2016_2022.csv`: coarse presence by year.
- `../1_raw_data/inputs/avg_winter_temp_yearly.csv`: Port Fourchon winter temperatures.

Outputs (`outputs/`)
- `filter1_species.csv`, `filter2_species.csv`, `filter3_species.csv`, `filter4_species.csv`, `filter_final.csv`: species retained after each filter stage.
- `plots/species_temp_plot.png`: CTI (weighted and unweighted) overlaid with winter temperature.

Software
- R >= 4.3 with tidyverse, ggplot2, ggrepel, janitor, brms, ggspatial, maptiles, sf.

Run order
1. Ensure STI/CTI outputs from Step 3 are available.
2. Run `ID_species_of_interest_BCODMO.R`.
3. Review filtered species lists and plot in `outputs/`.
