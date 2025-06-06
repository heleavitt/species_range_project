# Load required packages
library(rgbif)
library(tidyverse)
library(terra)
library(sp)
library(sdmpredictors)

setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")
options(
  gbif_user = Sys.getenv("GBIF_USER"),
  gbif_email = Sys.getenv("GBIF_EMAIL"),
  gbif_pwd = Sys.getenv("GBIF_PWD")
)


# Load presence data
presence_df <- read.csv("presence_pivot_merged_sp.csv")

# Calculate mean occurrence for each species
presence_df$mean_occurrence <- rowMeans(presence_df[, c("X2006", "X2016", "X2022_2023")], na.rm = TRUE)

# Filter to species-level taxa (space in name) or defined groups
species_level <- presence_df[grepl(" ", presence_df$Taxon), ]

# Unique taxa (including pooled groups)
target_taxa <- unique(species_level$Taxon)
# Define manual taxon groupings
group_taxa <- list(
  "Minuca spp." = c("Minuca longisignalis"),
  "Palaemon spp." = c("Palaemon pugio", "Palaemon vulgaris")
)

# Load SST raster
sst_layer <- load_layers("BO_sstmean", equalarea = FALSE)

# Convert raster and points to terra objects
sst_terra <- rast(sst_layer)  # convert from raster::raster to terra::SpatRaster

# Results container
results <- data.frame(Taxon = character(), STI = numeric(), IQR = numeric(), n_points = integer(), stringsAsFactors = FALSE)

write("", file = "gbif_citations.txt")  # clears or creates the file

# Loop through taxa
for (taxon in target_taxa) {
  message("Processing ", taxon)
  
  # Determine which species to query
  if (taxon %in% names(group_taxa)) {
    species_vec <- group_taxa[[taxon]]
  } else {
    species_vec <- taxon
  }
  for (sp in species_vec) {
    key <- name_backbone(name = sp)$usageKey
    
  # Submit the download request
  download_key <- occ_download(
    pred("taxonKey", key),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("basisOfRecord", "HUMAN_OBSERVATION"),
    pred_gte("year", 2005),
    format = "SIMPLE_CSV"
  )
  
  # Wait for GBIF to prepare the download
  message("Waiting for GBIF to finish download for ", sp, "...")
  
  # Wait loop
  while (occ_download_meta(download_key)$status != "SUCCEEDED") {
    Sys.sleep(2)  # wait 2 seconds
  }
  # Log DOI
  
  # Get metadata and citation
  meta <- occ_download_meta(download_key)
  
  # Log check
  print(meta)
  
  # catch issues with doi
  if (!is.null(meta$citation) && nzchar(meta$citation)) {
    cat(paste0(sp, ": ", meta$citation, "\n\n"), file = "gbif_citations.txt", append = TRUE)
  } else if (!is.null(meta$doi) && nzchar(meta$doi)) {
    cat(paste0(sp, ": https://doi.org/", meta$doi, "\n\n"), file = "gbif_citations.txt", append = TRUE)
  } else {
    warning("No DOI or citation found for ", sp)
    cat(paste0(sp, ": [No DOI found]\n\n"), file = "gbif_citations.txt", append = TRUE)
  }
    # Once ready, download and import
  if (!dir.exists("gbif_downloads")) dir.create("gbif_downloads")
  
  dl <- occ_download_get(download_key, 
                         path = "gbif_downloads", 
                         overwrite = TRUE)
  
  occ_data <- occ_download_import(dl)
  occ_all <- rbind(occ_data, occ_all)
  
  }
}

write.csv(occ_all, "occ_all.csv", row.names = FALSE)

  
  # Clean coordinates
  occ_clean <- occ_all %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
    distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  coords <- occ_clean %>% dplyr::select(decimalLongitude, decimalLatitude)
  sp_points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Extract SST with 5 km buffer
  pts <- vect(sp_points, crs = "EPSG:4326")
  
  # Extract SST using a 5 km buffer (terra uses degrees, so ~0.05 is 5km near equator)
  occ_clean$sst <- terra::extract(sst_terra, pts, buffer = 0.05, fun = mean)[, 2] 
  sst_vals <- na.omit(occ_clean$sst)

  species_STI <- occ_clean %>% drop_na(sst) %>% group_by(species) %>% summarise(sti = mean(sst),
                                                               sti_sd = sd(sst), .groups = 'keep')
  
  write.csv(species_STI, "STI_results_by_taxon.csv", row.names = FALSE)
  
  STI_Yearly<-presence_df %>% left_join(species_STI, by = c("Taxon" = "species")) %>% drop_na(sti) %>% mutate(sti2005 = sti*X2006,
                                                                                                              sti_2016 = sti*X2016,
                                                                                                              sti_2023 = sti*X2022_2023)
sum(STI_Yearly$sti2005)  
sum(STI_Yearly$sti_2016)  
sum(STI_Yearly$sti_2023)  
