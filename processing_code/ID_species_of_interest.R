library(tidyverse)
library(janitor)
library(brms)

setwd("C:/Users/hl51981/OneDrive - University of Georgia/Leavitt_Herbert/PFFW/Manuscripts/Global Change/Revision_repository")

pivot_clean<-read.csv("pivot_clean.csv", check.names = FALSE)
sti<-read.csv("STI_results_by_taxon.csv")

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
cor_threshold <- 0.5  # set your threshold
cor_results <- mean_abund_by_year %>%
  group_by(species) %>%
  summarise(
    n_years = n(),
    r = ifelse(n_years > 1, cor(Year, mean_abund, use = "complete.obs"), NA_real_),
    .groups = "drop"
  )
o6 <- cor_results %>% filter(r > .8 & r < 0.9)
o6_ploting <-abund_long %>% filter( species %in% o6$species) %>% group_by(Year, species) %>% summarise(abundance = mean(abundance))

ggplot(data = o6_ploting, aes(x = Year, y = log(abundance+1), colour = species))+
  geom_line() + 
  geom_point()
# Step 4: Determine presence by year (binary)
presence_matrix <- abund_long %>%
  mutate(present = ifelse(abundance > 0, 1, 0)) %>%
  group_by(species, Year) %>%
  summarise(present = max(present), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = present, values_fill = 0)

# Step 5: Apply logic filters
filter_flags <- presence_matrix %>%
  mutate(
    pass_presence = (`2016` == 1 & `2006` == 0) |
      (`2023` == 1 & `2006` == 0) |
      (`2023` == 1 & `2016` == 1 & `2006` == 0)
  )
# Step 6: Combine correlation + logical filter
pass_species <- cor_results %>%
  left_join(filter_flags, by = "species") %>%
  mutate(pass_filter1 = r > cor_threshold | pass_presence) %>%
  filter(pass_filter1)

# Output: species passing filter 1
co_migrator_candidates <- pass_species$species


##### Filter 2 #######
species_STI<-read.csv("STI_results_by_taxon.csv")
co_migrator_sti <- data.frame(species = co_migrator_candidates) %>% left_join(species_STI, by = "species")

co_migrator_filter2 <- co_migrator_sti %>% filter(sti > PtFou_vals)

# Filter 3 #### 
# for absences:
# Select relevant columns
target_data <- pivot_clean %>% 
  select(all_of(c(co_migrator_filter2$species, "Year", "SampleID")))

# Melt the dataset
long_data <- target_data %>%
  pivot_longer(cols = all_of(co_migrator_filter2$species), names_to = "species", values_to = "abundance") %>%
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
undetected_probs <- expand.grid(species = co_migrator_filter2$species, Year = unique(long_data$Year)) %>%
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

final_set<-detection_wide %>% filter(if_any(colnames(.[,c(5,7,8)]), ~. < .05))

# for increases 
consistent_species<-pass_species %>% filter(`2006` == 1 & `2016`== 1 & `2023` == 1)

for (species in consistent_species){
  a<-long_data %>% filter(species == species)
  a$Year_f <- factor(a$Year)
  mdl_a<-brm(
    abundance ~ Year_f,
    data = a,
    family = zero_inflated_negbinomial(),
    chains = 4, cores = 4, iter = 2500,
    control = list(adapt_delta = 0.995)
  )
  
  post<-as_draws_df(mdl_a)
  # 2016 vs baseline
  post$diff_2016_baseline <- post$b_Year_f2016
  
  # 2023 vs baseline
  post$diff_2023_baseline <- post$b_Year_f2023
  
  # 2023 vs 2016
  post$diff_2023_2016 <- post$b_Year_f2023 - post$b_Year_f2016
  
  print(paste(species, "likelihood of true increase from 2006 to 2016"))
  quantile(post$diff_2016_baseline, c(0.025, 0.5, 0.975))
  
  print(paste(species, "likelihood of true increase from 2006 to 2023"))
  quantile(post$diff_2023_baseline, c(0.025, 0.5, 0.975))
  
  print(paste(species, "likelihood of true increase from 2016 to 2023"))
  quantile(post$diff_2023_2016, c(0.025, 0.5, 0.975))
}
