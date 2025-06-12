library(readr)
library(rgbif)
library(tidyverse)
library(terra)
library(sp)
library(sdmpredictors)


setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")


# Set directory containing cleaned CSVs
clean_dir <- "gbif_downloads/clean_csvs"

# List all CSV files
csv_files <- list.files(clean_dir, pattern = "_clean\\.csv$", full.names = TRUE)

# Read and bind all files, adding species name as a column
gbif_list <- purrr::map_dfr(csv_files, function(file) {
  df <- read.csv(file)
  df$species_name <- gsub("_clean\\.csv$", "", basename(file)) %>%
    gsub("_", " ", .)
  return(df)
})

# Remove any NULLs from skipped folders
gbif_list <- Filter(Negate(is.null), gbif_list)

# Optional: Combine all data frames
combined_gbif <- do.call(rbind, gbif_list)


levels(as.factor(combined_gbif$species))
# Clean coordinates
occ_clean <- combined_gbif %>%
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
sp_points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract SST with 5 km buffer
pts <- vect(sp_points, crs = "EPSG:4326")

datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
layers <- list_layers(datasets)

# Load SST raster
sst_layer <- load_layers("BO_sstmin", equalarea = FALSE)

# Convert raster and points to terra objects
sst_terra <- rast(sst_layer)  # convert from raster::raster to terra::SpatRaster


# Extract SST using a 5 km buffer (terra uses degrees, so ~0.05 is 5km near equator)
occ_clean$sst <- terra::extract(sst_terra, pts, buffer = 0.05, fun = mean)[, 2] 
sst_vals <- na.omit(occ_clean$sst)


PtFou_Coords <- data.frame(y = 29.105560, x = -90.194443) 


PtFou_vec <- vect(PtFou_Coords, geom = c("x", "y"), crs = "EPSG:4326")

PtFou_vals <- terra::extract(sst_terra, PtFou_vec, buffer = 0.05, fun = mean)[, 2] 


# recompile split out species
species_STI <- occ_clean %>% drop_na(sst) %>% group_by(species) %>% summarise(sti = mean(sst),
                                                                              sti_sd = sd(sst), .groups = 'keep')

species_niche_size <- occ_clean %>% drop_na(sst) %>% drop_na(sal) %>% group_by(species) %>% summarise(temp_mean = mean(sst),
                                                                                                      temp_2.5 = quantile(sst, 0.025),
                                                                                                      temp_97.5 = quantile(sst, 0.975),
                                                                                                      temp_range = (temp_97.5-temp_2.5),
                                                                                                      .groups = "keep")
shrimpsst<-occ_clean %>% drop_na(sst) %>% drop_na(sal) %>%filter(species == "Penaeus setiferus") %>% select(sst)
hist( shrimpsst$sst)

write.csv(species_STI, "STI_results_by_taxon.csv", row.names = FALSE)

pivot_all<-read.csv("pivot_all.csv", check.names = FALSE)
sti_vector <- setNames(species_STI$sti, species_STI$species)

# Keep only species in both the matrix and STI vector
common_species <- intersect(names(sti_vector), colnames(pivot_all))

colnames(pivot_all)[! colnames(pivot_all) %in% (names(sti_vector))]
# Subset community matrix and STI vector to shared species
comm_sub <- as.matrix(pivot_all[, common_species])
rownames(comm_sub) <- pivot_all$SampleID  # or whatever holds your sample names

sti_sub <- sti_vector[common_species]

numerator <- comm_sub %*% sti_sub
denominator <- rowSums(comm_sub)
mean_sti <- as.numeric(numerator) / denominator

comm_df <- as.data.frame(comm_sub)  # convert back to data.frame
comm_df$SampleID <- rownames(comm_sub)  # add SampleID column
comm_df$mean_sti <- mean_sti  # attach calculated STI
comm_df$Year <- pivot_all$Year  # bring in Year from pivot_all

write.csv(comm_df, "pivot_clean.csv")



CTI_Yearly<-comm_df %>% drop_na(mean_sti) %>% group_by(Year) %>% summarize(mean_STI = mean(mean_sti),
                                                                           sd_STI = sd(mean_sti))
