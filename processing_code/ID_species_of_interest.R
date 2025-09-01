library(tidyverse)
library(janitor)
library(brms)
library(ggspatial)
library(ggrepel)
library(maptiles)
library(sf)


setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")

pivot_clean<-read.csv("outputs/pivot_clean.csv", check.names = FALSE)
sti<-read.csv("outputs/STI_results_by_taxon.csv")

pivot_clean <- pivot_clean[, names(pivot_clean) != ""]

pivot_clean <- pivot_clean %>%
  mutate(Year = ifelse(Year == "2022_2023", 2023, Year))

pivot_clean$Year <- as.numeric(pivot_clean$Year)

# Step 1: Reshape to long format
abund_long <- pivot_clean %>%
  pivot_longer(cols = -c(SampleID, Year), names_to = "species", values_to = "abundance")

# Step 2: Calculate mean abundance per year for each species
mean_abund_by_year <- abund_long %>%
  group_by(species, Year) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop")



# Step 3: Calculate correlation (abundance ~ year) for each species
cor_threshold <- 0.7  # set your threshold
cor_results <- mean_abund_by_year %>%
  group_by(species) %>%
  summarise(
    n_years = n(),
    r = ifelse(n_years > 1, cor(Year, mean_abund, use = "complete.obs"), NA_real_),
    .groups = "drop"
  )


# Step 4: Determine presence by year (binary)
presence_matrix <- abund_long %>%
  mutate(present = ifelse(abundance > 0, 1, 0)) %>%
  group_by(species, Year) %>%
  summarise(present = max(present), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = present, values_fill = 0)

# Step 5: Apply logic filters
filter_flags <- presence_matrix %>%
  mutate(
    pass_presence =       (`2023` == 1 & `2006` == 0) |
      (`2023` == 1 & `2016` == 1 & `2006` == 0)
  )
# Step 6: Combine correlation + logical filter
pass_species <- cor_results %>%
  left_join(filter_flags, by = "species") %>%
  mutate(pass_filter1 = r > cor_threshold | pass_presence) %>%
  filter(pass_filter1)

# Output: species passing filter 1: Shows an increase in abundance from 2006 - 2023. 
co_migrator_candidates <- pass_species$species
write.csv(co_migrator_candidates, "outputs/filter1_species.csv")
##### Filter 2: Species has enough high-quality gbif observationsw/in grind cells to confidently draw conclusions (n > 5 grids)
species_STI<-read.csv("outputs/STI_results_by_taxon.csv")
species_STI <- species_STI %>% filter(n_kept > 5)
co_migrator_sti <- data.frame(species = co_migrator_candidates) %>% left_join(species_STI, by = "species") %>% na.omit()
write.csv(co_migrator_sti, "outputs/filter2_species.csv")

##### Filter 3: Species that have an Minimum temperature (actually defined as the 2.5th percintile of bioclim temps) 
##### greater than 3C, which is approcimately the minimum winter temp of Port FOurchon in 2004 #######

co_migrator_filter3 <- co_migrator_sti %>% filter(sti_2.5 > 3) %>% left_join(cor_results, by = "species")
write.csv(co_migrator_filter3, "outputs/filter3_species.csv")

##### Filter 4 Species that have are near the northern-most extend of their observed range(<1% of observations north of 31 degrees)
co_migrator_filter4 <- co_migrator_filter3 %>% filter(northern_lat < 30.7) 
write.csv(co_migrator_filter4, "outputs/filter4_species.csv")
# Filter 5 #### 
# Calculates the likelyhood of false absences based on the commonness of the species and the number of samples taken each year. 
# returns species with a less than 5% of false absence in an indevidual year or 
# over several years. 

# Select relevant columns
target_data <- pivot_clean %>% 
  select(all_of(c(co_migrator_filter4$species, "Year", "SampleID")))

# Melt the dataset
long_data <- target_data %>%
  pivot_longer(cols = all_of(co_migrator_filter4$species), names_to = "species", values_to = "abundance") %>%
  mutate(detected = abundance > 0)

samples_per_year <- long_data %>%
  distinct(Year, SampleID) %>%
  count(Year, name = "sites")

# Estimate average detection probability for each species (when detected)
detection_probs <- long_data %>%
  filter(detected) %>%
  group_by(species, Year) %>%
  summarise(samples_detected = n(), .groups = 'drop') %>%
  right_join(samples_per_year, by = "Year") %>%
  group_by(species) %>% na.omit() %>% 
  summarise(detection_prob = mean(samples_detected / sites), .groups = 'drop')


# Merge and calculate P(Undetected | Present)
undetected_probs <- expand.grid(species = co_migrator_filter4$species, Year = unique(long_data$Year)) %>%
  left_join(long_data %>% group_by(species, Year) %>% summarise(detected = any(detected), .groups = 'drop'),
            by = c("species", "Year")) %>%
  left_join(detection_probs, by = "species") %>%
  left_join(samples_per_year, by = "Year") %>%
  mutate(p_false_absence = ifelse(!detected, (1 - detection_prob)^sites, NA))

# Pivot to wide format
detection_wide <- undetected_probs %>%
  pivot_wider(
    id_cols = species,
    names_from = Year,
    values_from = c(detected, p_false_absence),
    names_glue = "{Year}_{.value}"
  )

# Add joint probability for 2006 and 2016 if undetected in both
detection_wide <- detection_wide %>%
  mutate(p_false_absence_2006_2016 = ifelse(
    `2006_detected` == FALSE & `2016_detected` == FALSE,
    `2006_p_false_absence` * `2016_p_false_absence`,
    NA
  ))

# View results
print(detection_wide)


# filter 
### THIS TABEL REPRESENTS THE SET OF SPECIES PASSING ALL FILTERS, 
###### FIlter 1 1: 
possible_set<-detection_wide %>% filter(if_any(colnames(.[,c(5,7,8)]), ~. < .05))

write.csv(possible_set, "outputs/filter_final.csv")



# for increases library(dplyr)
library(brms)
library(posterior)
library(tibble)

co_migrator_filter2<-co_migrator_filter3 %>% left_join(pass_species, by = "species")

# Filter species present in all three years
consistent_species <- co_migrator_filter2 %>%
  filter(`2006` == 1 & `2016` == 1 & `2023` == 1) %>%
  pull(species)

# Create empty list to collect results
summary_list <- list()

for (sp in consistent_species) {
  message("Running model for: ", sp)
  
  a <- long_data %>% filter(species == sp)
  a$Year_f <- factor(a$Year)
  
  # Fit model
  mdl <- brm(
    abundance ~ Year_f,
    data = a,
    family = zero_inflated_negbinomial(),
    chains = 4, cores = 4, iter = 2000,
    control = list(adapt_delta = 0.995),
    silent = TRUE
  )
  
  # Extract posterior draws
  post <- as_draws_df(mdl)
  
  # Compute differences
  post$diff_2016_baseline <- post$b_Year_f2016
  post$diff_2023_baseline <- post$b_Year_f2023
  post$diff_2023_2016 <- post$b_Year_f2023 - post$b_Year_f2016
  
  # Summarize into credible intervals and medians
  summary_row <- tibble(
    species = sp,
    
    # 2016 vs 2006
    diff_2016_q025 = quantile(post$diff_2016_baseline, 0.025),
    diff_2016_median = quantile(post$diff_2016_baseline, 0.5),
    diff_2016_q975 = quantile(post$diff_2016_baseline, 0.975),
    
    # 2023 vs 2006
    diff_2023_q025 = quantile(post$diff_2023_baseline, 0.025),
    diff_2023_median = quantile(post$diff_2023_baseline, 0.5),
    diff_2023_q975 = quantile(post$diff_2023_baseline, 0.975),
    
    # 2023 vs 2016
    diff_23_16_q025 = quantile(post$diff_2023_2016, 0.025),
    diff_23_16_median = quantile(post$diff_2023_2016, 0.5),
    diff_23_16_q975 = quantile(post$diff_2023_2016, 0.975)
  )
  
  # Store
  summary_list[[sp]] <- summary_row
}

# Combine all summaries into a single data frame
species_model_summary <- bind_rows(summary_list)

possible_set2 <-  species_model_summary %>%
  filter(if_all(-species, ~ . > 0)) %>%
  pull(species)

final_set <- c(possible_set$species, possible_set2)

# Save to CSV if needed
write.csv(species_model_summary, "species_model_summary.csv", row.names = FALSE)


######### BINOMIAL REGRESSION 
presence <- read.csv("outputs/species_presence_2006_2016_2022.csv")

presence_long <- presence %>%
  pivot_longer(cols = c("X2005", "X2016", "X2022"),
               names_to = "Year", values_to = "present") %>%
  filter(present == TRUE)

presence_long$Year <- as.numeric(gsub("X", "", presence_long$Year))



# Join with STI data using valid_name
presence_sti <- presence_long %>%
  left_join(sti, by = c("valid_name" = "species")) %>%
  filter(!is.na(sti)) 

# Calculate CTI by year
presence_cti <- presence_sti %>%
  group_by(Year) %>%
  summarise(CTI = mean(sti_2.5, na.rm = TRUE))%>%  # remove unmatched species
  mutate(Type = "presence_cti")

cti_fauna_yearly <- pivot_clean %>%
  drop_na(mean_sti) %>%
  group_by(Year) %>%
  summarize(fauna_cti = mean(mean_sti),
            fauna_min_cti = mean(mean_sti_2.5),
            sd_fauna_min = sd(mean_sti_2.5),
            sd_fauna_cti = sd(mean_sti), 
            fauna_range = mean(mean_sti_range)) %>%
  mutate(Year = gsub("2022_2023", "2022", Year))


cti_fauna_yearly$Year <- as.numeric(cti_fauna_yearly$Year)

mangrove_sti <- species_STI[species_STI$species == "Avicennia germinans", "sti_2.5"]
marsh_sti <- species_STI[species_STI$species == "Sporobolus alterniflorus", "sti_2.5"]


cti_flora_yearly <- data.frame(Year = c(2002, 2014, 2022),
                               mangrove_ptl = c(.0231, .0961, .1891),
                               marsh_ptl = c(0.7947, 0.7333, 0.610)) %>%
  mutate(veg_cti = (mangrove_ptl * mangrove_sti) + (marsh_ptl*marsh_sti)
  )

########### PLOTTING ########################
# --- ensure climate columns are named consistently ---
climate <- read.csv("outputs/avg_winter_temp_yearly.csv") %>%
  dplyr::rename(Year = year, temp = X0)

# --- build CTI (exclude vegetation) and recode labels for legend ---
combined_cti <- dplyr::bind_rows(
  cti_fauna_yearly %>%
    dplyr::select(Year, CTI = fauna_min_cti) %>%
    dplyr::mutate(Type = "Weighted\nCommunity\nTemperature\nIndex"),
  presence_cti %>%
    dplyr::transmute(Year, CTI, Type = "Unweighted\nCommunity\nTemperature\nIndex"), 
  cti_flora_yearly %>%  mutate( Type = "Average Min.\nTemp of\nVegetation\nCover") %>% dplyr::select(Year, CTI = veg_cti, Type)
)

combined_cti[4,"Year"] <- 2006
combined_cti[3, "Year"] <- 2022

# --- compute a pure offset so climate overlays CTI on the primary axis ---
# Align by overall means across overlapping years
yrs <- intersect(combined_cti$Year, climate$Year)
offset <- mean(combined_cti$CTI[combined_cti$Year %in% yrs], na.rm = TRUE) -
  mean(climate$temp[climate$Year %in% yrs], na.rm = TRUE)

climate_plot <- climate %>%
  dplyr::mutate(temp_shifted = temp + offset)

# palette & shapes
col_vals <- c(
  "Weighted\nCommunity\nTemperature\nIndex" = "#ee3377",
  "Average Min.\nTemp of\nVegetation\nCover" = "forestgreen",
  "Unweighted\nCommunity\nTemperature\nIndex"     = "#EE7733",
  "Port Fourchon\nTemperature"        = "black"
)
shape_vals <- c(
  "Weighted Community\nTemperature Index" = 16,  # circle
  "Unweighted Community\ntemperature Index"     = 18   # diamond
)

library(ggplot2)
library(ggrepel)



png("outputs/plots/species_temp_plot.png", width = 18, height = 20, units = "cm", res = 600)
ggplot() +
  # CTI
  geom_line(data = combined_cti,
            aes(Year, CTI, color = Type), linewidth = 1.5) +
  geom_point(data = combined_cti,
             aes(Year, CTI, color = Type), size = 3.5) +
  # Climate (color only)
  geom_smooth(
    data = climate_plot,
    aes(Year, temp_shifted, color = "Port Fourchon\nTemperature"),
    method = "loess", size = 1.5, se = TRUE, linetype = "dashed"
  ) +
  
  # --- one shared legend title ---
  scale_color_manual(
    name   = "",
    values = col_vals,
    breaks = c("Unweighted\nCommunity\nTemperature\nIndex",
               "Weighted\nCommunity\nTemperature\nIndex",
               "Average Min.\nTemp of\nVegetation\nCover",
               "Port Fourchon\nTemperature")
  )  +
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1)) +
  # --------------------------------

scale_y_continuous(
  name = "Community Temperature Index (°C)",
  sec.axis = sec_axis(~ . - offset,
                      name = "Port Fourchon winter temperature (°C)")
) +
  scale_x_continuous(breaks = c(2006, 2015, 2022)) +
  labs(x = "Year", title = "Figure 3") +
  theme_bw() +
  guides(       # <-- its own layer
    color = guide_legend(order = 1, nrow = 1),
    shape = guide_legend(order = 1, nrow = 1)
  ) +
  theme(
    legend.position = "bottom",        # move legend to bottom
    legend.key.width = unit(1.5, "cm"),
    legend.text  = element_text(size = 10),
    axis.title.y.right = element_text(margin = margin(l = 10))
  )

dev.off()
######### extra plotting code for including species if you wanted #########

# species STI reference lines (on CTI axis)
geom_segment(data = filter(species_STI, species %in% final_set |
                             species %in% c("Penaeus setiferus",
                                            "Avicennia germinans",
                                            "Fundulus grandis",
                                            "Callinectes sapidus",
                                            "Armases americanum",
                                            "Palaemon spp.",
                                            "Sporobolus alterniflorus")),
             aes(x = 2002, xend = 2023, y = sti_2.5, yend = sti_2.5),
             linetype = "dotted", color = "gray40") +
  
  geom_text_repel(
    data = filter(species_STI, species %in% final_set |
                    species %in% c("Penaeus setiferus",
                                   "Avicennia germinans",
                                   "Fundulus grandis",
                                   "Callinectes sapidus",
                                   "Armases americanum",
                                   "Palaemon spp.",
                                   "Sporobolus alterniflorus")),
    aes(x = 2023, y = sti_2.5, label = species),
    nudge_x = 3, direction = "y", hjust = 0, size = 3,
    fontface = "italic", color = "gray40",
    segment.color = "gray70", segment.size = 0.3
  ) +









# Recode trend as binary
species_STI <- species_STI %>% filter(trend != "no change") %>% 
  mutate(trend_bin = ifelse(trend == "negative", 1, 0)) %>%
  drop_na(trend_bin, sti)

# Run binomial regression
trend_model <- glm(trend_bin ~ sti + sti_range, data = species_STI_r, family = binomial())
# View summary
summary(trend_model)
plot(trend_model)

ggplot(species_STI_r, aes(x = sti, y = trend_bin)) +
  geom_jitter(width = 0.1, height = 0.05, alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  theme_minimal() +
  labs(x = "Species Temperature Index (STI)",
       y = "Probability of Decline",
       title = "STI vs. Abundance Trend (Binary)")

# View summary
summary(trend_model)
plot(trend_model)

options(
  gbif_user = Sys.getenv("GBIF_USER"),
  gbif_email = Sys.getenv("GBIF_EMAIL"),
  gbif_pwd = Sys.getenv("GBIF_PWD")
)

  for (sp in final_set) {
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
    dl <- occ_download_get(download_key, path = "gbif_downloads/final_set", overwrite = TRUE)
    
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
      filter(coordinateUncertaintyInMeters < 10000 | is.na(coordinateUncertaintyInMeters)) %>%
      filter(!coordinateUncertaintyInMeters %in% c(301, 3036, 999, 9999)) %>%
      filter(!(decimalLatitude == 0 | decimalLongitude == 0)) %>%
      cc_cen(buffer = 2000) %>%
      cc_cap(buffer = 2000) %>%
      cc_inst(buffer = 2000) %>%
      distinct(decimalLongitude, decimalLatitude, speciesKey, datasetKey, .keep_all = TRUE)
    
    # Save cleaned data
    csv_path <- paste0("gbif_downloads/final_set/final_set_csvs/", gsub(" ", "_", sp), "_final.csv")
    write.csv(clean_data, csv_path, row.names = FALSE)
    
    # Log citation or DOI
    meta <- occ_download_meta(download_key)
    if (!is.null(meta$citation) && nzchar(meta$citation)) {
      cat(paste0(sp, ": ", meta$citation, "\n\n"), file = "gbif_final_set_citations.txt", append = TRUE)
    } else if (!is.null(meta$doi) && nzchar(meta$doi)) {
      cat(paste0(sp, ": https://doi.org/", meta$doi, "\n\n"), file = "gbif_final_set_citations.txt", append = TRUE)
    } else {
      warning("No DOI or citation found for ", sp)
      cat(paste0(sp, ": [No DOI found]\n\n"), file = "gbif_final_set_citations.txt", append = TRUE)
    }
  }






# Set directory containing cleaned CSVs
clean_dir <- "gbif_downloads/final_set/final_set_csvs/"

# List all CSV files
csv_files <- list.files(clean_dir, pattern = "_final\\.csv$", full.names = TRUE)

final_set_list <- purrr::map_dfr(csv_files, function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  df$catalogNumber <- as.character(df$catalogNumber)
  df$recordNumber <- as.character(df$recordNumber)
  df$species_name <- gsub("_final.csv$", "", basename(file)) %>%
    gsub("_", " ", .)
  return(df)
})


north_obs_sf <- st_as_sf(final_set_list, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

north_obs_sf_filtered <- north_obs_sf %>%
  filter(st_coordinates(.)[, 2] > 20)

# Download basemap tile (choose provider: e.g., "CartoDB.Positron", "OpenStreetMap")
basemap <- get_tiles(north_obs_sf_filtered, provider = "CartoDB.Positron", crop = TRUE)

# Plot with ggplot2
ggplot() +
  layer_spatial(basemap) +
  geom_sf(data = north_obs_sf_filtered, aes(color = species), size = 1.2, alpha = 0.2) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl") +
  theme_minimal() +
  labs(title = "GBIF Observations North of Port Fourchon (1980–2005)",
       color = "Species")# Optionally: write to CSV for inspection
write.csv(north_obs_sf_filtered, "GBIF_North_Observations_1980_2005.csv", row.names = FALSE)

# Extract species flagged as false positives
false_positives <- unique(north_obs_df$species)

# Filter final set
final_set_cleaned <- setdiff(final_set, false_positives)
