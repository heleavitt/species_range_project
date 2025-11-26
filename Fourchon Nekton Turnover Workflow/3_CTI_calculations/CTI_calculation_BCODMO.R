# Community Temperature Index (CTI) and Species Thermal Index (STI)
# Annotated script for BCODMO submission
#
# Purpose:
#   - Read cleaned GBIF occurrence records for focal taxa.
#   - Extract mean sea surface temperature (SST) for each occurrence.
#   - Thin spatially dense records to reduce sampling bias.
#   - Derive per-species thermal metrics (STI and confidence bounds).
#   - Combine with community matrices to compute CTI-style site metrics.
#
# Key inputs (relative to repository root):
#   - 2_gbif_workflow/clean_csvs/*_clean.csv  : cleaned GBIF occurrences per taxon.
#   - raw_data/species_presence_2006_2016_2022.csv : presence metadata used in joins.
#   - outputs/pivot_all.csv                 : community matrix with SampleID + Year.
#
# Key outputs:
#   - outputs/species_latitudes.csv         : latitudinal summaries for each taxon.
#   - outputs/STI_results_by_taxon.csv      : STI statistics per taxon after thinning.
#   - outputs/pivot_clean.csv               : community matrix with CTI metrics appended.
#
# Notes:
#   - Script assumes WGS84 coordinates and uses a ~5 km buffer when extracting SST.
#   - Spatial thinning uses 50 km hexagonal grid (adjust `cell_km` if needed).
#   - All file paths are set relative to the repository; update `setwd()` if moving.

library(readr)          # CSV I/O helpers
library(rgbif)          # GBIF data structures (cleaned data already supplied)
library(tidyverse)      # Data wrangling (dplyr, purrr, etc.)
library(terra)          # Raster handling for SST extraction
library(sf)             # Modern simple features for spatial thinning
library(sp)             # Legacy spatial classes (used for initial point creation)
library(ggrepel)        # Loaded for plotting if needed (not used directly here)
library(sdmpredictors)  # Access to environmental rasters (WorldClimSST)

# Set working directory to repository root so relative paths resolve.
setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")

# --------------------------------------------------------------------------- #
# 1. Load cleaned GBIF occurrence data
# --------------------------------------------------------------------------- #

clean_dir <- "gbif_downloads/clean_csvs"
csv_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)

# Bind all per-species CSVs into one data frame, adding a `species_name` column.
gbif_list <- purrr::map_dfr(csv_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df$catalogNumber <- as.character(df$catalogNumber)
  df$recordNumber  <- as.character(df$recordNumber)
  df$species_name  <- gsub("_clean\\.csv$", "", basename(file)) %>% gsub("_", " ", .)
  return(df)
})

# --------------------------------------------------------------------------- #
# 2. Clean coordinates, normalize taxon naming, and prepare spatial objects
# --------------------------------------------------------------------------- #

occ_clean <- gbif_list %>%
  # Remove rows without coordinates and deduplicate identical coordinate pairs
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
  # Cast numeric columns
  mutate(across(
    c(decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters,
      coordinatePrecision, elevation, elevationAccuracy, depth, depthAccuracy),
    ~ as.numeric(.)
  )) %>%
  # Collapse selected taxa into genus-level buckets for analysis
  # Minuca was not identified to species in 3/4 seasons, but longisignalis is by far the most prevalent at our sites, so we are simplifying all Minuca observations to be longisignalis 
  # Palaemon was not identified to species in one of the seasons. Pugio and vulgaris are both common and have similar ranges, so they were compiled for this analysis
  # Callinecres, Alpheus, and Farfantepenaus were often identified to genus level in 2005 and 2006. Subsiquent years found that these genuses were overwhelmingly attibuted to one species (> 90%)
  # so these genus-level observations were attributed to the most common species in each genus (C. sapidus, A. heterochaelis, Penaus (Farfantepenaus) aztecus)
  # Panopeus herbstii was a species complex that was split into Panopeus obesus and Panopeus simpsoni, the new convention was used in 2022 and 2023 but not other years, so P. simpsoni and P. obesus were compiled 
  # into Panopeus spp. and analyzed together. These species have similar ranges, but the modern P. herbstii is mostly confined to the atlantic coast. 
  # Species will be compiled to their "analysis ready" form here so that statistics can be run on processing outputs. 
  mutate(
    species = case_when(
      species == "Minuca longisignalis" ~ "Minuca spp.",
      species == "Palaemon pugio"       ~ "Palaemon spp.",
      species == "Palaemon vulgaris"    ~ "Palaemon spp.",
      species == "Panopeus simpsoni"    ~ "Panopeus spp.",
      species == "Panopeus obesus"      ~ "Panopeus spp.",
      species == "Panopeus herbstii"    ~ "Panopeus spp.",
      TRUE                              ~ species
    )
  ) %>%
  as.data.frame()

# Build SpatialPoints for SST extraction (WGS84)
coords <- occ_clean %>%
  dplyr::select(decimalLongitude, decimalLatitude) %>%
  as.matrix()
sp_points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))

# --------------------------------------------------------------------------- #
# 3. Retrieve SST layer and extract mean SST per occurrence (5 km buffer)
# --------------------------------------------------------------------------- #

datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
layers   <- list_layers(datasets)

# WorldClim bioclimatic variable 6 is mean temperature of coldest month (proxy SST)
sst_layer <- load_layers("WC_bio6", equalarea = FALSE)
layer_citations(layers = c("WC_bio6"), astext = FALSE)  # Keep citation on hand

# Convert to terra objects and extract buffered means (~0.05 deg ≈ 5 km at equator)
sst_terra <- rast(sst_layer)
pts       <- vect(sp_points, crs = "EPSG:4326")

occ_clean$sst <- terra::extract(sst_terra, pts, buffer = 0.05, fun = mean)[, 2]
sst_vals      <- na.omit(occ_clean$sst)  # retained for quick QA if needed

# Reference point: Port Fourchon mean SST (same buffer)
PtFou_Coords <- data.frame(y = 29.105560, x = -90.194443)
PtFou_vec    <- vect(PtFou_Coords, geom = c("x", "y"), crs = "EPSG:4326")
PtFou_vals   <- terra::extract(sst_terra, PtFou_vec, buffer = 0.05, fun = mean)[, 2]

# --------------------------------------------------------------------------- #
# 4. Spatial thinning to reduce sampling bias
# --------------------------------------------------------------------------- #

cell_km <- 50  # Grid size (km) for thinning; increase to thin more aggressively.

# Keep only usable rows for thinning
occ_use <- occ_clean %>%
  drop_na(sst, decimalLongitude, decimalLatitude, species)

# Convert to sf, project to meters (EPSG:3857) for distance-based operations
pts_sf <- st_as_sf(occ_use,
                   coords = c("decimalLongitude", "decimalLatitude"),
                   crs = 4326, remove = FALSE)
pts_m  <- st_transform(pts_sf, 3857)

# Create a hexagonal grid covering all points
grid_hex <- st_make_grid(st_as_sfc(st_bbox(pts_m)),
                         cellsize = cell_km * 1000,
                         square = FALSE)
grid_hex <- st_sf(grid_id = seq_along(grid_hex), geometry = grid_hex)

# Assign each point to a hex cell; if a point falls on a boundary, attach nearest
pts_thin <- st_join(pts_m, grid_hex, left = TRUE, join = st_within)
idx_na <- which(is.na(pts_thin$grid_id))
if (length(idx_na)) {
  nearest <- st_nearest_feature(pts_m[idx_na, ], grid_hex)
  pts_thin$grid_id[idx_na] <- grid_hex$grid_id[nearest]
}

# Keep one random record per species per grid cell
thin <- pts_thin |>
  st_drop_geometry() |>
  group_by(species, grid_id) |>
  slice_sample(n = 1) |>
  ungroup()

# --------------------------------------------------------------------------- #
# 5. Latitudinal summaries (per-species min/mean/max)
# --------------------------------------------------------------------------- #

sp_presence <- read.csv("./raw_data/species_presence_2006_2016_2022.csv")

species_latitudes <- thin %>%
  filter(decimalLongitude >= -130 & decimalLongitude <= -30) %>%
  group_by(species) %>%
  summarize(
    "Min Latitude"  = min(decimalLatitude),
    "Mean Latitude" = mean(decimalLatitude),
    "Max Latitude"  = max(decimalLatitude)
  ) %>%
  left_join(sp_presence, by = c("species" = "valid_name")) %>%
  na.omit()

write.csv(species_latitudes, "outputs/species_latitudes.csv", row.names = FALSE)

# --------------------------------------------------------------------------- #
# 6. Species Thermal Index (STI) metrics after thinning
# --------------------------------------------------------------------------- #

species_STI <- thin %>%
  group_by(species) %>%
  summarise(
    sti          = mean(sst, na.rm = TRUE),
    sti_2.5      = quantile(sst, 0.025, na.rm = TRUE),
    sti_97.5     = quantile(sst, 0.975, na.rm = TRUE),
    sti_range    = sti_97.5 - sti_2.5,
    northern_lat = quantile(decimalLatitude, 0.99, na.rm = TRUE),  # robust max
    n_kept       = dplyr::n(),
    .groups      = "drop"
  )

write.csv(species_STI, "outputs/STI_results_by_taxon.csv", row.names = FALSE)

# --------------------------------------------------------------------------- #
# 7. Community Temperature Index (CTI) calculations
# --------------------------------------------------------------------------- #

pivot_all <- read.csv("outputs/pivot_all.csv", check.names = FALSE)

# Prepare lookup vectors keyed by species name
sti_vector       <- setNames(species_STI$sti, species_STI$species)
sti_2.5_vector   <- setNames(species_STI$sti_2.5, species_STI$species)
sti_range_vector <- setNames(species_STI$sti_range, species_STI$species)

# Align species between community matrix and STI results
common_species <- intersect(names(sti_vector), colnames(pivot_all))
comm_sub       <- as.matrix(pivot_all[, common_species])
rownames(comm_sub) <- pivot_all$SampleID

sti_sub       <- sti_vector[common_species]
sti_2.5_sub   <- sti_2.5_vector[common_species]
sti_range_sub <- sti_range_vector[common_species]

# Weighted means: matrix multiplication of community counts by STI vectors
numerator        <- comm_sub %*% sti_sub
numerator_2.5    <- comm_sub %*% sti_2.5_sub
numerator_range  <- comm_sub %*% sti_range_sub
denominator      <- rowSums(comm_sub)

mean_sti        <- as.numeric(numerator) / denominator
mean_sti_2.5    <- as.numeric(numerator_2.5) / denominator
mean_sti_range  <- as.numeric(numerator_range / denominator)

# Reattach to community data frame with metadata
comm_df <- as.data.frame(comm_sub)
comm_df$SampleID       <- rownames(comm_sub)
comm_df$mean_sti       <- mean_sti
comm_df$mean_sti_2.5   <- mean_sti_2.5
comm_df$Year           <- pivot_all$Year
comm_df$mean_sti_range <- mean_sti_range

write.csv(comm_df, "outputs/pivot_clean.csv", row.names = FALSE)

# End of script
