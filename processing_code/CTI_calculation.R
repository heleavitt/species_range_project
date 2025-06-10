library(readr)

zip_files <- list.files("gbif_downloads", pattern = "\\.zip$", full.names = TRUE)

unzipped_base <- "gbif_downloads/unzipped_data"

# Create a folder to extract to
dir.create(unzipped_base, showWarnings = FALSE)

# Unzip all
lapply(zip_files, function(zf) {
  unzip(zf, exdir = file.path(unzipped_base, tools::file_path_sans_ext(basename(zf))))
})

# Define base folder containing all unzipped subfolders

# Get a list of all subfolders
subfolders <- list.dirs(unzipped_base, full.names = TRUE, recursive = FALSE)

# For each subfolder, find the first .csv file and read it

gbif_list <- lapply(subfolders, function(subdir) {
  csv_file <- list.files(subdir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_file) == 0) return(NULL)
  message("Reading: ", csv_file[1])
  tryCatch(
    read_delim(csv_file[1], delim = "\t",col_types = cols(.default = "c") ,show_col_types = FALSE),
    error = function(e) {
      warning("Skipping ", csv_file[1], ": ", e$message)
      return(NULL)
    }
  )
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
      
      TRUE ~ species  # leave all others unchanged
    )
  )


coords <- occ_clean %>%
  dplyr::select(decimalLongitude, decimalLatitude) %>%
  as.matrix()
sp_points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))

# Extract SST with 5 km buffer
pts <- vect(sp_points, crs = "EPSG:4326")

# Extract SST using a 5 km buffer (terra uses degrees, so ~0.05 is 5km near equator)
occ_clean$sst <- terra::extract(sst_terra, pts, buffer = 0.05, fun = mean)[, 2] 
sst_vals <- na.omit(occ_clean$sst)

# recompile split out species
species_STI <- occ_clean %>% drop_na(sst) %>% group_by(species) %>% summarise(sti = mean(sst),
                                                                              sti_sd = sd(sst), .groups = 'keep')

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





STI_Yearly<-comm_df %>% drop_na(mean_sti) %>% group_by(Year) %>% summarize(mean_STI = mean(mean_sti),
                                                                              sd_STI = sd(mean_sti))
