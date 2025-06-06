# Load required packages
library(rgbif)
library(dplyr)
library(raster)
library(sp)
library(sdmpredictors)

options(
  gbif_user = Sys.getenv("GBIF_USER"),
  gbif_email = Sys.getenv("GBIF_EMAIL"),
  gbif_pwd = Sys.getenv("GBIF_PWD")
)


# Load presence data
presence_df <- read.csv("presence_pivot.csv")

# Calculate mean occurrence for each species
presence_df$mean_occurrence <- rowMeans(presence_df[, c("X2006", "X2016", "X2022_2023")], na.rm = TRUE)

# Filter to species-level taxa (space in name) or defined groups
species_level <- presence_df[grepl(" ", presence_df$Taxon), ]

# Rename species into pooled groupings
species_level$Taxon[species_level$Taxon %in% c("Minuca", "Minuca longisignalis", "Minuca pugnax", "Minuca rapax")] <- "Minuca spp."
species_level$Taxon[species_level$Taxon %in% c("Leander tenuicornis", "Palaemon intermedius", "Palaemon pugio", "Palaemon vulgaris", "Palaemonetes")] <- "Palaemon spp."

# Unique taxa (including pooled groups)
target_taxa <- unique(species_level$Taxon)
# Define manual taxon groupings
group_taxa <- list(
  "Minuca spp." = c("Minuca longisignalis"),
  "Palaemon spp." = c("Palaemon pugio", "Palaemon vulgaris")
)

# Load SST raster
sst_layer <- load_layers("BO_sstmean", equalarea = FALSE)

# Results container
results <- data.frame(Taxon = character(), STI = numeric(), IQR = numeric(), n_points = integer(), stringsAsFactors = FALSE)

# Loop through taxa
for (taxon in target_taxa) {
  message("Processing ", taxon)
  
  # Determine which species to query
  if (taxon %in% names(group_taxa)) {
    species_vec <- group_taxa[[taxon]]
  } else {
    species_vec <- taxon
  }
  
  # Get occurrence data for all species in group
  occ_all <- data.frame()
  for (sp in species_vec) {
    key <- name_backbone(name = sp)$usageKey
    download_key <- occ_download(
      pred("taxonKey", key),
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      pred("basisOfRecord", "HUMAN_OBSERVATION"),
      pred_gte("year", 2005),
      format = "SIMPLE_CSV"
    )
    
    # Wait until ready, or query later:
    message("Download started for ", sp, ": ", download_key)
    Sys.sleep(60)  # add manual delay or loop back later
    
    dl <- occ_download_get(download_key, overwrite = TRUE)
    occ_data <- occ_download_import(dl)
    occ_all <- bind_rows(occ_all, occ_data)
  }
  
  # Clean coordinates
  occ_clean <- occ_all %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
    distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  coords <- occ_clean %>% select(decimalLongitude, decimalLatitude)
  sp_points <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Extract SST with 5 km buffer
  occ_clean$sst <- extract(sst_layer, sp_points, buffer = 5000, fun = mean, na.rm = TRUE)
  sst_vals <- na.omit(occ_clean$sst)
  
  # Skip if not enough data
  if (length(sst_vals) < 10) {
    warning("Too few points for ", taxon)
    next
  }
  
  # Compute stats
  sti <- mean(sst_vals)
  iqr <- IQR(sst_vals)
  
  # Save results
  results <- rbind(results, data.frame(Taxon = taxon, STI = sti, IQR = iqr, n_points = length(sst_vals)))
}

# Output
print(results)
write.csv(results, "STI_results_by_taxon.csv", row.names = FALSE)