library(readr)
library(rgbif)
library(tidyverse)
library(terra)
library(sf)
library(sp)
library(ggrepel)

library(sdmpredictors)


setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")


# Set directory containing cleaned CSVs
clean_dir <- "gbif_downloads/clean_csvs"

# List all CSV files
csv_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)

gbif_list <- purrr::map_dfr(csv_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df$catalogNumber <- as.character(df$catalogNumber)
  df$recordNumber <- as.character(df$recordNumber)
  df$species_name <- gsub("_clean\\.csv$", "", basename(file)) %>%
    gsub("_", " ", .)
  return(df)
})

unique(gbif_list$species)
# Clean coordinates
occ_clean <- gbif_list %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>% 
  mutate(across(
    c(decimalLatitude, decimalLongitude, coordinateUncertaintyInMeters,
      coordinatePrecision, elevation, elevationAccuracy,
      depth, depthAccuracy),
    ~ as.numeric(.)
  )) %>% as.data.frame()%>%
  mutate( # recomple split out species 
    species = case_when(
      species == "Minuca longisignalis" ~ "Minuca spp.",
      species == "Palaemon pugio" ~ "Palaemon spp.",
      species == "Palaemon vulgaris" ~ "Palaemon spp.",
      species == "Panopeus simpsoni" ~ "Panopeus spp.",
      species == "Panopeus obesus" ~ "Panopeus spp.",
      species == "Panopeus herbstii" ~ "Panopeus spp.",

      TRUE ~ species  # leave all others unchanged
    )
  )

coords <- occ_clean %>%
  dplyr::select(decimalLongitude, decimalLatitude) %>%
  as.matrix()

sp_points <- SpatialPoints(coords,
                           proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract SST with 5 km buffer
pts <- vect(sp_points, crs = "EPSG:4326")

datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
layers <- list_layers(datasets)

# Load SST raster
sst_layer <- load_layers("WC_bio6", equalarea = FALSE)
layer_citations(layers = c("WC_bio6"), astext = FALSE)
# Convert raster and points to terra objects
sst_terra <- rast(sst_layer)  # convert from raster::raster to terra::SpatRaster


# Extract SST using a 5 km buffer 
# (terra uses degrees, so ~0.05 is 5km near equator)
occ_clean$sst <- terra::extract(sst_terra, pts, buffer = 0.05, fun = mean)[, 2] 
sst_vals <- na.omit(occ_clean$sst)


PtFou_Coords <- data.frame(y = 29.105560, x = -90.194443) 


PtFou_vec <- vect(PtFou_Coords, geom = c("x", "y"), crs = "EPSG:4326")

PtFou_vals <- terra::extract(sst_terra,
                             PtFou_vec, buffer = 0.05, fun = mean)[, 2]


########## Spatial thinning to reduce bias from high number of samples in certain
########## regions
cell_km <- 50   # ≈ thinning grid size (km). Increase to thin more aggressively.

# Keep usable rows
occ_use <- occ_clean %>%
  drop_na(sst, decimalLongitude, decimalLatitude, species)

# To sf; project to meters for distance-based thinning
pts <- st_as_sf(occ_use,
                coords = c("decimalLongitude","decimalLatitude"),
                crs = 4326, remove = FALSE)
pts_m <- st_transform(pts, 3857)

# Hex grid & IDs
grid_hex <- st_make_grid(st_as_sfc(st_bbox(pts_m)),
                         cellsize = cell_km * 1000,
                         square = FALSE)
grid_hex <- st_sf(grid_id = seq_along(grid_hex), geometry = grid_hex)

# Assign each point to a hex cell (use within; repair edge NAs with nearest)
pts_thin <- st_join(pts_m, grid_hex, left = TRUE, join = st_within)
idx_na <- which(is.na(pts_thin$grid_id))
if (length(idx_na)) {
  nearest <- st_nearest_feature(pts_m[idx_na, ], grid_hex)
  pts_thin$grid_id[idx_na] <- grid_hex$grid_id[nearest]
}

# ≤1 record per species per cell
thin <- pts_thin |>
  st_drop_geometry() |>
  group_by(species, grid_id) |>
  slice_sample(n = 1) |>
  ungroup()


sp_presence<-read.csv("./raw_data/species_presence_2006_2016_2022.csv")

species_latitudes<-thin %>%
  filter(decimalLongitude >= -130 & decimalLongitude <= -30) %>% 
  group_by(species) %>% summarize("Min Latitude" = min(decimalLatitude), 
                                  "Mean Latitude" = mean(decimalLatitude),
                                  "Max Latitude" = max(decimalLatitude)) %>% 
  left_join(sp_presence, by = c( "species"= "valid_name")) %>% na.omit()
write.csv(species_latitudes, "outputs/species_latitudes.csv")

# Recompute STI metrics on the thinned set
species_STI <- thin %>%
  group_by(species) %>%
  summarise(
    sti         = mean(sst, na.rm = TRUE),
    sti_2.5     = quantile(sst, 0.025, na.rm = TRUE),
    sti_97.5    = quantile(sst, 0.975, na.rm = TRUE),
    sti_range   = sti_97.5 - sti_2.5,
    northern_lat = quantile(decimalLatitude, 0.99, na.rm = TRUE), # robust-ish max
    n_kept      = dplyr::n(),
    .groups     = "drop"
  )

write.csv(species_STI, "outputs/STI_results_by_taxon.csv", row.names = FALSE)

pivot_all<-read.csv("outputs/pivot_all.csv", check.names = FALSE)
sti_vector <- setNames(species_STI$sti, species_STI$species)
sti_2.5_vector <-setNames(species_STI$sti_2.5, species_STI$species)
sti_range_vector <-setNames(species_STI$sti_range, species_STI$species)

# Keep only species in both the matrix and STI vector
common_species <- intersect(names(sti_vector), colnames(pivot_all))

colnames(pivot_all)[! colnames(pivot_all) %in% (names(sti_vector))]
# Subset community matrix and STI vector to shared species
comm_sub <- as.matrix(pivot_all[, common_species])
rownames(comm_sub) <- pivot_all$SampleID  # or whatever holds your sample names

sti_sub <- sti_vector[common_species]
sti_2.5_sub <- sti_2.5_vector[common_species]
sti_range_sub <- sti_range_vector[common_species]



numerator <- comm_sub %*% sti_sub
numerator_2.5 <- comm_sub %*% sti_2.5_sub
numerator_range <- comm_sub %*% sti_range_sub

denominator <- rowSums(comm_sub)
mean_sti <- as.numeric(numerator) / denominator
mean_sti_2.5 <- as.numeric(numerator_2.5) / denominator
mean_sti_range <- as.numeric(numerator_range  / denominator)


comm_df <- as.data.frame(comm_sub)  # convert back to data.frame
comm_df$SampleID <- rownames(comm_sub)  # add SampleID column
comm_df$mean_sti <- mean_sti  # attach calculated STI
comm_df$mean_sti_2.5 <-mean_sti_2.5
comm_df$Year <- pivot_all$Year  # bring in Year from pivot_all
comm_df$mean_sti_range <- mean_sti_range
write.csv(comm_df, "outputs/pivot_clean.csv")






