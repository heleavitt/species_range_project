# GBIF Download & Cleaning (BCODMO annotated)
#
# Purpose:
#   - Read list of target taxa from presence summaries.
#   - Download GBIF occurrence data (post-2005, human observations with coordinates).
#   - Apply coordinate/quality filters and basic spatial cleaning (CoordinateCleaner).
#   - Save cleaned per-taxon CSVs and log GBIF citations/DOIs.
#
# Key inputs (relative to repo root):
#   - outputs/presence_pivot_merged_sp.csv : list of taxa to query.
#   - Environment variables for GBIF API: GBIF_USER, GBIF_EMAIL, GBIF_PWD.
#
# Key outputs:
#   - gbif_downloads/clean_csvs/<taxon>_clean.csv : cleaned occurrences per taxon.
#   - gbif_citations.txt                          : citation/DOI log for downloads.
#
# Notes:
#   - Threatened taxa use a larger coordinate uncertainty ceiling (30 km).
#   - Grouped taxa (e.g., Minuca spp.) are expanded to species for queries.
#   - CoordinateCleaner flags (cap, cen, inst) applied with 2 km buffer.

library(rgbif)
library(tidyverse)
library(terra)
library(sp)
library(sdmpredictors)
library(CoordinateCleaner)

# Set working directory to repository root
setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")

# Configure GBIF credentials from env vars
options(
  gbif_user  = Sys.getenv("GBIF_USER"),
  gbif_email = Sys.getenv("GBIF_EMAIL"),
  gbif_pwd   = Sys.getenv("GBIF_PWD")
)

# --------------------------------------------------------------------------- #
# 1) Build target taxon list
# --------------------------------------------------------------------------- #

presence_df <- read.csv("outputs/presence_pivot_merged_sp.csv")

# Species-level taxa have a space in the name; keep pooled groups explicitly
species_level <- presence_df %>% filter(str_detect(Taxon, " "))

target_taxa <- unique(species_level$Taxon) %>%
  c("Avicennia germinans", "Spartina alterniflora")

# Expand pooled taxa into species queries
group_taxa <- list(
  "Minuca spp."   = c("Minuca longisignalis"),
  "Palaemon spp." = c("Palaemon pugio", "Palaemon vulgaris"),
  "Panopeus spp." = c("Panopeus obesus", "Panopeus simpsoni")
)

# Use wider uncertainty allowance for threatened taxa
threatened_taxa <- c("Fundulus jenkinsi", "Negaprion brevirostris")

# Ensure output directory exists
dir.create("gbif_downloads/clean_csvs", showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------- #
# 2) Download, clean, and save per-taxon GBIF data
# --------------------------------------------------------------------------- #

for (taxon in target_taxa) {
  message("Processing ", taxon)

  # Larger uncertainty threshold for threatened taxa
  uncert <- if (taxon %in% threatened_taxa) 30000 else 1000

  # Determine which species to query (expand pooled taxa)
  species_vec <- if (taxon %in% names(group_taxa)) group_taxa[[taxon]] else taxon

  for (sp in species_vec) {
    key <- name_backbone(name = sp)$usageKey

    # Submit download: human observations, coordinates present, post-2005
    download_key <- occ_download(
      pred("taxonKey", key),
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      pred("basisOfRecord", "HUMAN_OBSERVATION"),
      pred_gte("year", 2005),
      format = "SIMPLE_CSV"
    )

    occ_download_wait(download_key)

    # Fetch and import
    dl <- occ_download_get(download_key, path = "gbif_downloads", overwrite = TRUE)
    clean_data <- occ_download_import(dl)

    # Skip if nothing returned or no coordinates
    if (nrow(clean_data) == 0 || !"decimalLatitude" %in% names(clean_data)) {
      warning("No usable data returned for ", sp)
      next
    }

    # Core cleaning steps
    clean_data <- clean_data %>%
      filter(occurrenceStatus == "PRESENT") %>%
      filter(!basisOfRecord %in% c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN")) %>%
      filter(year >= 1900) %>%
      filter(coordinatePrecision < 0.01 | is.na(coordinatePrecision)) %>%
      filter(coordinateUncertaintyInMeters < uncert | is.na(coordinateUncertaintyInMeters)) %>%
      filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 9999)) %>%
      filter(!(decimalLatitude == 0 | decimalLongitude == 0)) %>%
      cc_cen(buffer = 2000) %>%   # centroid check
      cc_cap(buffer = 2000) %>%   # capitals
      cc_inst(buffer = 2000) %>%  # institutions
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

# End of script
