# Species of Interest Identification (BCODMO annotated)
#
# Purpose:
#   - Use CTI outputs to identify co-migrator candidates showing increasing abundance.
#   - Apply sequential filters: abundance trend, GBIF data sufficiency, thermal limits,
#     range position, and false-absence likelihood.
#   - Produce summary CTI plots against winter temperature.
#
# Key inputs (relative to repo root):
#   - 3_CTI_calculations\outputs\pivot_clean.csv                                                  : community matrix with CTI metrics attached.
#   - 3_CTI_calculations\outputs\STI_results_by_taxon.csv                                         : species thermal indices from GBIF extractions.
#   - 1_raw_data\inputs\species_presence_2006_2016_2022.csv/species_presence_2006_2016_2022.csv   : coarse presence by year.
#   - raw_data/avg_winter_temp_yearly.csv                                                         : winter temperature time series (Port Fourchon).
#
# Key outputs:
#   - outputs/filter1_species.csv
#   - outputs/filter2_species.csv
#   - outputs/filter3_species.csv
#   - outputs/filter4_species.csv
#   - outputs/filter_final.csv            : species passing all filters.
#   - outputs/plots/species_temp_plot.png : CTI vs winter temperature figure.

library(tidyverse)
library(janitor)
library(brms)
library(ggspatial)
library(ggrepel)
library(maptiles)
library(sf)
library(ggplot2)


# --------------------------------------------------------------------------- #
# 1) Load inputs
# --------------------------------------------------------------------------- #

pivot_clean <- read.csv("3_CTI_calculations/outputs/pivot_clean.csv", check.names = FALSE)
sti         <- read.csv("3_CTI_calculations/outputs/STI_results_by_taxon.csv")

# Remove empty column names if present
pivot_clean <- pivot_clean[, names(pivot_clean) != ""]

# Normalize Year coding (treat 2022/2023 combined as 2023)
pivot_clean <- pivot_clean %>%
  mutate(Year = ifelse(Year == "2022_2023", 2023, Year)) %>%
  mutate(Year = as.numeric(Year))

# --------------------------------------------------------------------------- #
# 2) Abundance trend and presence logic (Filters 1–4)
# --------------------------------------------------------------------------- #

# Long format for abundance across years
abund_long <- pivot_clean %>%
  pivot_longer(cols = -c(SampleID, Year), names_to = "species", values_to = "abundance")

# Mean abundance per year per species
mean_abund_by_year <- abund_long %>%
  group_by(species, Year) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop")

# Correlation of abundance vs year (Filter 1: positive trend or presence gain)
cor_threshold <- 0.7
cor_results <- mean_abund_by_year %>%
  group_by(species) %>%
  summarise(
    n_years = n(),
    r = ifelse(n_years > 1, cor(Year, mean_abund, use = "complete.obs"), NA_real_),
    .groups = "drop"
  )

# Presence matrix by year (binary)
presence_matrix <- abund_long %>%
  mutate(present = ifelse(abundance > 0, 1, 0)) %>%
  group_by(species, Year) %>%
  summarise(present = max(present), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = present, values_fill = 0)

# Logical presence gain criteria (detected recently, absent earlier)
filter_flags <- presence_matrix %>%
  mutate(
    pass_presence = (`2023` == 1 & `2006` == 0) |
                    (`2023` == 1 & `2016` == 1 & `2006` == 0)
  )

# Combine correlation and presence gains
pass_species <- cor_results %>%
  left_join(filter_flags, by = "species") %>%
  mutate(pass_filter1 = r > cor_threshold | pass_presence) %>%
  filter(pass_filter1)

co_migrator_candidates <- pass_species$species
write.csv(co_migrator_candidates, "outputs/filter1_species.csv", row.names = FALSE)

# Filter 2: require sufficient GBIF grid cells (n_kept > 5)
species_STI <- read.csv("3_CTI_calculations/outputs/STI_results_by_taxon.csv")
species_STI <- species_STI %>% filter(n_kept > 5)

co_migrator_sti <- data.frame(species = co_migrator_candidates) %>%
  left_join(species_STI, by = "species") %>%
  na.omit()
write.csv(co_migrator_sti, "outputs/filter2_species.csv", row.names = FALSE)

# Filter 3: minimum thermal tolerance above ~3 C (Port Fourchon winter baseline)
co_migrator_filter3 <- co_migrator_sti %>%
  filter(sti_2.5 > 3) %>%
  left_join(cor_results, by = "species")
write.csv(co_migrator_filter3, "outputs/filter3_species.csv", row.names = FALSE)

# Filter 4: near northern edge of observed range (99th percentile lat < 30.7)
co_migrator_filter4 <- co_migrator_filter3 %>%
  filter(northern_lat < 30.7)
write.csv(co_migrator_filter4, "outputs/filter4_species.csv", row.names = FALSE)

# --------------------------------------------------------------------------- #
# 3) False-absence probability (Filter 5)
# --------------------------------------------------------------------------- #

# Subset community data to candidate species
target_data <- pivot_clean %>%
  select(all_of(c(co_migrator_filter4$species, "Year", "SampleID")))

# Long format with detection flags
long_data <- target_data %>%
  pivot_longer(cols = all_of(co_migrator_filter4$species),
               names_to = "species", values_to = "abundance") %>%
  mutate(detected = abundance > 0)

# Sites sampled per year
samples_per_year <- long_data %>%
  distinct(Year, SampleID) %>%
  count(Year, name = "sites")

# Average detection probability per species (when detected)
detection_probs <- long_data %>%
  filter(detected) %>%
  group_by(species, Year) %>%
  summarise(samples_detected = n(), .groups = "drop") %>%
  right_join(samples_per_year, by = "Year") %>%
  group_by(species) %>%
  na.omit() %>%
  summarise(detection_prob = mean(samples_detected / sites), .groups = "drop")

# Probability of false absence by year
undetected_probs <- expand.grid(species = co_migrator_filter4$species,
                                Year = unique(long_data$Year)) %>%
  left_join(long_data %>%
              group_by(species, Year) %>%
              summarise(detected = any(detected), .groups = "drop"),
            by = c("species", "Year")) %>%
  left_join(detection_probs, by = "species") %>%
  left_join(samples_per_year, by = "Year") %>%
  mutate(p_false_absence = ifelse(!detected, (1 - detection_prob)^sites, NA))

# Wide table summarizing detection and false-absence probabilities
detection_wide <- undetected_probs %>%
  pivot_wider(
    id_cols = species,
    names_from = Year,
    values_from = c(detected, p_false_absence),
    names_glue = "{Year}_{.value}"
  ) %>%
  mutate(
    p_false_absence_2006_2016 = ifelse(
      `2006_detected` == FALSE & `2016_detected` == FALSE,
      `2006_p_false_absence` * `2016_p_false_absence`,
      NA
    )
  )

print(detection_wide)

# Filter: keep species with <5% false-absence probability in any key comparison
possible_set <- detection_wide %>%
  filter(if_any(colnames(.[, c(5, 7, 8)]), ~ . < 0.05))

write.csv(possible_set, "outputs/filter_final.csv", row.names = FALSE)

# --------------------------------------------------------------------------- #
# 4) CTI calculations and plot against winter temperature
# --------------------------------------------------------------------------- #

presence <- read.csv("1_raw_data/inputs/species_presence_2006_2016_2022.csv")

presence_long <- presence %>%
  pivot_longer(cols = c("X2005", "X2016", "X2022"),
               names_to = "Year", values_to = "present") %>%
  filter(present == TRUE) %>%
  mutate(Year = as.numeric(gsub("X", "", Year)))

presence_sti <- presence_long %>%
  left_join(sti, by = c("valid_name" = "species")) %>%
  filter(!is.na(sti))

presence_cti <- presence_sti %>%
  group_by(Year) %>%
  summarise(CTI = mean(sti_2.5, na.rm = TRUE), .groups = "drop") %>%
  mutate(Type = "Unweighted\nCommunity\nTemperature\nIndex")

cti_fauna_yearly <- pivot_clean %>%
  drop_na(mean_sti) %>%
  group_by(Year) %>%
  summarize(
    fauna_cti      = mean(mean_sti),
    fauna_min_cti  = mean(mean_sti_2.5),
    sd_fauna_min   = sd(mean_sti_2.5),
    sd_fauna_cti   = sd(mean_sti),
    fauna_range    = mean(mean_sti_range),
    .groups = "drop"
  ) %>%
  mutate(Year = gsub("2022_2023", "2022", Year)) %>%
  mutate(Year = as.numeric(Year))

# Flora thermal contributions (fixed proportions by year)
mangrove_sti <- species_STI[species_STI$species == "Avicennia germinans", "sti_2.5"]
marsh_sti    <- species_STI[species_STI$species == "Sporobolus alterniflorus", "sti_2.5"]

cti_flora_yearly <- data.frame(
  Year       = c(2002, 2014, 2022),
  mangrove_ptl = c(0.0231, 0.0961, 0.1891),
  marsh_ptl    = c(0.7947, 0.7333, 0.610)
) %>%
  mutate(veg_cti = (mangrove_ptl * mangrove_sti) + (marsh_ptl * marsh_sti)) %>%
  mutate(Type = "Average Min.\nTemp of\nVegetation\nCover") %>%
  select(Year, CTI = veg_cti, Type)

# Combine CTI series
combined_cti <- bind_rows(
  cti_fauna_yearly %>%
    select(Year, CTI = fauna_min_cti) %>%
    mutate(Type = "Weighted\nCommunity\nTemperature\nIndex"),
  presence_cti,
  cti_flora_yearly
)

# Ensure years align for plotting labels
combined_cti[combined_cti$Type == "Weighted\nCommunity\nTemperature\nIndex" & combined_cti$Year == 2022, "Year"] <- 2022
combined_cti[combined_cti$Type == "Weighted\nCommunity\nTemperature\nIndex" & combined_cti$Year == 2006, "Year"] <- 2006

# Climate time series and offset to align scales
climate <- read.csv("raw_data/avg_winter_temp_yearly.csv") %>%
  rename(Year = year, temp = X0)

yrs <- intersect(combined_cti$Year, climate$Year)
offset <- mean(combined_cti$CTI[combined_cti$Year %in% yrs], na.rm = TRUE) -
  mean(climate$temp[climate$Year %in% yrs], na.rm = TRUE)

climate_plot <- climate %>%
  mutate(temp_shifted = temp + offset)

# Colors and shapes
col_vals <- c(
  "Weighted\nCommunity\nTemperature\nIndex" = "#ee3377",
  "Average Min.\nTemp of\nVegetation\nCover" = "forestgreen",
  "Unweighted\nCommunity\nTemperature\nIndex" = "#EE7733",
  "Port Fourchon\nTemperature" = "black"
)

png("outputs/plots/species_temp_plot.png", width = 18, height = 20, units = "cm", res = 600)
ggplot() +
  geom_line(data = combined_cti, aes(Year, CTI, color = Type), linewidth = 1.5) +
  geom_point(data = combined_cti, aes(Year, CTI, color = Type), size = 3.5) +
  geom_smooth(
    data = climate_plot,
    aes(Year, temp_shifted, color = "Port Fourchon\nTemperature"),
    method = "loess", size = 1.5, se = TRUE, linetype = "dashed"
  ) +
  scale_color_manual(
    name   = "",
    values = col_vals,
    breaks = c("Unweighted\nCommunity\nTemperature\nIndex",
               "Weighted\nCommunity\nTemperature\nIndex",
               "Average Min.\nTemp of\nVegetation\nCover",
               "Port Fourchon\nTemperature")
  ) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 1)) +
  scale_y_continuous(
    name = "Community Temperature Index (deg C)",
    sec.axis = sec_axis(~ . - offset, name = "Port Fourchon winter temperature (deg C)")
  ) +
  scale_x_continuous(breaks = c(2006, 2015, 2022)) +
  labs(x = "Year", title = "Figure 3") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 10),
    axis.title.y.right = element_text(margin = margin(l = 10))
  )
dev.off()

# End of script
