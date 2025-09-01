library(vegan)
library(tidyverse)


pivot_clean<-read.csv("pivot_clean.csv", check.names = FALSE)

pivot_clean <- pivot_clean[, names(pivot_clean) != ""]

pivot_clean <- pivot_clean %>%
  mutate(Year = ifelse(Year == "2022_2023", 2023, Year))

# Extract site and year info
metadata <- pivot_clean %>% select(SampleID, Year, mean_sti, mean_sti_2.5, mean_sti_range)

# Extract only species data (assumes 3rd column onward are species)
comm_matrix <- pivot_clean %>% select(-SampleID, -Year, -mean_sti, -mean_sti_2.5, -mean_sti_range)


# Optional: check for zeros or highly sparse data
summary(rowSums(comm_matrix))  # Should not be all zero

# Calculate row sums
row_totals <- rowSums(comm_matrix)

# Keep only samples with non-zero totals
valid_rows <- row_totals > 0
comm_matrix_clean <- comm_matrix[valid_rows, ]
metadata_clean <- metadata[valid_rows, ]

comm_hellinger <- decostand(comm_matrix, method = "hellinger")


nmds_result <- metaMDS(comm_hellinger, distance = "bray", k = 2, trymax = 500)
nmds_coords <- as.data.frame(scores(nmds_result, display = "sites"))

cti_fit <- envfit(nmds_result, metadata %>% select(mean_sti, mean_sti_range), permutations = 999)
plot(nmds_result, display = "sites")
plot(cti_fit, p.max = 0.05, col = "blue")

print(cti_fit)

# Force alignment
metadata_clean <- metadata_clean[rownames(site_scores), ]
site_scores$SampleID <- rownames(site_scores)
site_scores <- left_join(site_scores, metadata_clean, by = "SampleID")
metadata$SampleID <- as.character(metadata$SampleID)

nmds_coords <- left_join(nmds_coords, metadata, by = "SampleID")

top_nmds1 <- species_scores_df %>%
  arrange(desc(abs(NMDS1))) %>%
  slice_head(n = 10)
print(top_nmds1)

##### envfit
# Extract arrow vector for CTI
cti_arrow <- as.data.frame(cti_fit$vectors$arrows)
cti_arrow$Variable <- rownames(cti_arrow)

# Optionally scale arrow for better display
arrow_scale <- 1.5
cti_arrow <- cti_arrow %>%
  mutate(xend = NMDS1 * arrow_scale,
         yend = NMDS2 * arrow_scale)

site_scores <- as.data.frame(scores(nmds_result, display = "sites"))
site_scores$SampleID <- metadata$SampleID

# Join with metadata
site_scores <- left_join(site_scores, metadata, by = "SampleID")



ggplot(site_scores, aes(x = NMDS1, y = NMDS2, color = as.factor(Year))) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Year), linetype = 2) +
  # Add CTI arrow
  geom_segment(data = cti_arrow,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm")),
               inherit.aes = FALSE,
               color = "black") +
  geom_text(data = cti_arrow,
            aes(x = xend * 1.1, y = yend * 1.1, label = Variable),
            inherit.aes = FALSE,
            size = 4, fontface = "bold") +
  theme_minimal() +
  labs(title = "NMDS of Community Composition with CTI Vector",
       color = "Year")


dispersion <- betadisper(vegdist(comm_hellinger), metadata_clean$Year)
anova(dispersion)  # Confirmed significant
library(permute)
library(pairwiseAdonis)

# Extract distance matrix and group variable
distance_matrix <- vegdist(comm_hellinger)
groups <- metadata_clean$Year

# Run pairwise tests for homogeneity of dispersion (multivariate spread)
pairwise_result <- pairwise.perm.manova(distance_matrix, groups, nperm = 999)
print(pairwise_result)


adonis_cti <- adonis2(comm_hellinger ~ Year, data = metadata_clean, method = "bray", permutations = 999)
print(adonis_cti)
