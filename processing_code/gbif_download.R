# Load required packages
library(rgbif)
library(tidyverse)
library(terra)
library(sp)
library(sdmpredictors)
library(CoordinateCleaner)

setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")
options(
  gbif_user = Sys.getenv("GBIF_USER"),
  gbif_email = Sys.getenv("GBIF_EMAIL"),
  gbif_pwd = Sys.getenv("GBIF_PWD")
)


# Load presence data
presence_df<- read.csv("outputs/presence_pivot_merged_sp.csv")

# Filter to species-level taxa (space in name) or defined groups
species_level <- presence_df_new[grepl(" ", presence_df_new$Taxon), ]

# Unique taxa (including pooled groups)
target_taxa <- unique(species_level$Taxon) %>% rbind("Avicennia germinans", "Spartina alterniflora")


group_taxa <- list(
  "Minuca spp." = c("Minuca longisignalis"),
  "Palaemon spp." = c("Palaemon pugio", "Palaemon vulgaris"),
  "Panopeus spp." = c("Panopeus obesus", "Panopeus simpsoni")
)

threatened_taxa <- c("Fundulus jenkinsi", "Negaprion brevirostris") 

# Loop through taxa
for (taxon in target_taxa) {
  message("Processing ", taxon)
  if(taxon %in% threatened_taxa){
    uncert <- 30000
  }else{uncert <- 1000}
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
    
    occ_download_wait(download_key) 
    
    # Save raw zip
    dl <- occ_download_get(download_key, path = "gbif_downloads", overwrite = TRUE)
    
    clean_data <- occ_download_import(dl)
    

    # Skip species if no records are returned
    if (nrow(clean_data) == 0) {
      warning("No data returned for ", sp)
      next
    }
  
    if (nrow(clean_data) == 0 || !"decimalLatitude" %in% names(clean_data)) {
      warning("No usable data returned for ", sp)
      next
    }
    
    
    # Clean the data
    clean_data <- clean_data %>%
      filter(occurrenceStatus == "PRESENT") %>%
      filter(!basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")) %>%
      filter(year >= 1900) %>%
      filter(coordinatePrecision < 0.01 | is.na(coordinatePrecision)) %>%
      filter(coordinateUncertaintyInMeters < 30000 | is.na(coordinateUncertaintyInMeters)) %>%
      filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 9999)) %>%
      filter(!(decimalLatitude == 0 | decimalLongitude == 0)) %>%
      cc_cen(buffer = 2000) %>%
      cc_cap(buffer = 2000) %>%
      cc_inst(buffer = 2000) %>%
      distinct(decimalLongitude, decimalLatitude, speciesKey, datasetKey, .keep_all = TRUE)
    
    # Save cleaned data
    csv_path <- paste0("gbif_downloads/clean_csvs/", gsub(" ", "_", sp), "_clean.csv")
    write.csv(clean_data, csv_path, row.names = FALSE)
    
    # Log citation or DOI
    meta <- occ_download_meta(download_key)
    if (!is.null(meta$citation) && nzchar(meta$citation)) {
      cat(paste0(sp, ": ", meta$citation, "\n\n"), file = "gbif_citations.txt", append = TRUE)
    } else if (!is.null(meta$doi) && nzchar(meta$doi)) {
      cat(paste0(sp, ": https://doi.org/", meta$doi, "\n\n"), file = "gbif_citations.txt", append = TRUE)
    } else {
      warning("No DOI or citation found for ", sp)
      cat(paste0(sp, ": [No DOI found]\n\n"), file = "gbif_citations.txt", append = TRUE)
    }
  }
}

