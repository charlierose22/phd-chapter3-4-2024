
# LOCATION ----------------------------------------------------------------

# Create a list of target antibiotics.
target_antibiotics <- unique(assay_samples$target_antibiotics_major)

# split based on target antibiotics for location 
split_location <- split(means_loc, means_loc$target_antibiotics_major)
loc_amino <- split_location$Aminoglycoside
loc_beta <- split_location$`Beta Lactam`
loc_int <- split_location$Integron
loc_mdr <- split_location$MDR
loc_mlsb <- split_location$MLSB
loc_mge <- split_location$MGE
loc_other <- split_location$Other
loc_phen <- split_location$Phenicol
loc_sulf <- split_location$Sulfonamide
loc_tet <- split_location$Tetracycline
loc_quin <- split_location$Quinolone
loc_vanc <- split_location$Vancomycin
loc_trim <- split_location$Trimethoprim

# split based on target antibiotics for time 
split_time <- split(means_time, means_time$target_antibiotics_major)
time_amino <- split_time$Aminoglycoside
time_beta <- split_time$`Beta Lactam`
time_int <- split_time$Integron
time_mdr <- split_time$MDR
time_mlsb <- split_time$MLSB
time_mge <- split_time$MGE
time_other <- split_time$Other
time_phen <- split_time$Phenicol
time_sulf <- split_time$Sulfonamide
time_tet <- split_time$Tetracycline
time_quin <- split_time$Quinolone
time_vanc <- split_time$Vancomycin
time_trim <- split_time$Trimethoprim

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data <- location_study[location_study$target_antibiotics_major == 
                           target_antibiotic, ]
  
  # Create a heatmap of the data.
  ggplot(data, aes(x = day, y = gene, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "gene", colour = "abundance") +
    facet_grid(length ~ height) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("arg-analysis/figures/heatmaps/location-heatmap-", 
                target_antibiotic, ".png"), width = 7, height = 7)
}

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data2 <- location_study[location_study$target_antibiotics_major == 
                            target_antibiotic, ]
  
  # Create bar graphs of the data.
  data2 %>% 
    group_by(length, height) %>% 
    ggplot(aes(x = day, y = mean, fill = gene)) +
    geom_bar(position = "stack", stat = "identity") +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = gene), width = 0.2) +
    scale_fill_viridis(discrete = T) +
    scale_color_viridis(discrete = T) +
    labs(x = "day", y = "abundance", fill = "gene") +
    facet_grid(length ~ height) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("arg-analysis/figures/bargraph/location-bargraph-", 
                target_antibiotic, ".png"), width = 7, height = 7)
}

# TIME --------------------------------------------------------------------

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data3 <- time_study[time_study$target_antibiotics_major == 
                        target_antibiotic, ]
  
  # Create a heatmap of the data.
  ggplot(data3, aes(x = day, y = gene, fill = mean)) +
    geom_tile() +
    scale_y_discrete(limits = rev) +
    scale_fill_viridis(discrete = F) +
    labs(x = "day", y = "gene", colour = "abundance") +
    theme_ipsum(base_size = 10)
  
  # Save the heatmap to a file.
  ggsave(paste0("arg-analysis/figures/heatmaps/time-heatmap-", 
                target_antibiotic, ".png"), width = 6, height = 5)
}

# Create a loop for each target antibiotic.
for (target_antibiotic in target_antibiotics) {
  
  # Create a subset of the data for the current target antibiotic.
  data4 <- time_study[time_study$target_antibiotics_major == 
                        target_antibiotic, ]
  
  # Create line graphs of the data.
  ggplot(data4, aes(x = day, 
                    y = mean, 
                    color = gene)) +
    geom_point() +
    geom_errorbar(aes(x = day,
                      ymin = mean - se,
                      ymax = mean + se),
                  width = .6) +
    geom_line(aes(color =  gene)) +
    labs(x = "day", y = "abundance", color = "gene") +
    scale_color_viridis(discrete = T) +
    theme_ipsum(base_size = 10)
  
  # Save the linegraph to a file.
  ggsave(paste0("arg-analysis/figures/linegraph/time-error-linegraph-", 
                target_antibiotic, ".png"), width = 6, height = 5)
}

# 16S ---------------------------------------------------------------------

# create a graph for total abundance of 16S data across samples.
ggplot(means_16S, aes(x = day, y = mean, fill = height)) +
  geom_col(position = 'dodge') +
  geom_errorbar(aes(x = day,
                    y = mean,
                    ymin = mean - sd, 
                    ymax = mean + sd), 
                width = 1) +
  labs(x = "day", y = "ct", fill = "location") +
  scale_color_viridis(discrete = F) +
  theme_ipsum(base_size = 10)
ggsave(paste0("arg-analysis/figures/16S-november.png"), width = 8, height = 5)


# INDIVIDUAL --------------------------------------------------------------
library(wesanderson)
# aminoglycoside Location
loc_amino %>%
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_viridis(discrete = T) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_ipsum(base_size = 10)

# beta-lactam Location
loc_beta %>%
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  scale_fill_manual(values = wes_palette("GrandBudapest2", n = 2)) +
  theme_ipsum(base_size = 12)

loc_amino
loc_int
loc_mdr
loc_mlsb
loc_mge
loc_other
loc_phen
loc_sulf
loc_tet
loc_quin
loc_vanc
loc_trim

library(RColorBrewer)

time_beta %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "abundance", color = "gene") +
  ylim(0, 0.00015) +
  scale_color_manual(values = brewer.pal("Set3", n = 10)) +
  theme_ipsum(base_size = 12)
