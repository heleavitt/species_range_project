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


######### BINOMIAL REGRESSION 
presence <- read.csv("raw_data/species_presence_2006_2016_2022.csv")

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
climate <- read.csv("raw_data/avg_winter_temp_yearly.csv") %>%
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
