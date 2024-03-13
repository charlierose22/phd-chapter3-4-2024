
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
library(RColorBrewer)
# aminoglycoside Location
loc_amino %>%
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

# beta-lactam Location
loc_beta %>%
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  theme_minimal(base_size = 12)

loc_int %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_mdr %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_mlsb %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_mge %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_other %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_phen %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_sulf %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_tet %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_quin %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_vanc %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)

loc_trim %>% 
  ggplot(aes(x = height, y = mean, fill = day)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.6)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2,
                position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = brewer.pal("Dark2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "day") +
  facet_wrap(~gene, scales = "free") +
  theme_minimal(base_size = 12)


library(RColorBrewer)
darkpalette <- c("mediumblue", 
                 "darkred", 
                 "darkorange4", 
                 "orange", 
                 "darkgreen", 
                 "green", 
                 "darkmagenta", 
                 "hotpink", 
                 "darkturquoise", 
                 "tomato")
# with cpha
time_beta %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)
# zoomed in on bottom for no cpha
time_beta %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  ylim(0, 0.00015) +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_amino %>% 
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_int %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_mdr %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_mlsb %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_mge %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_other %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_phen %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_sulf %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_tet %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_quin %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_vanc %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

time_trim %>%
  ggplot(aes(x = day, y = mean)) +
  geom_point(aes(color = gene)) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se,
                    color = gene),
                width = .6) +
  geom_line(aes(color =  gene)) +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)


# -----
means_loc2 <- means_loc

means_loc2$std = NULL
means_loc2$se = NULL
means_loc2$n = NULL

loc_means_wide <- pivot_wider(means_loc2, names_from = day, values_from = mean)
colnames(loc_means_wide) <- c("gene", "height", "class", "day_1", "day_29")
loc_means_wide <- loc_means_wide[!grepl('MGE|Other|Integrons|MDR', loc_means_wide$class),]

loc_means_wide %>%
  group_by(class) %>%
  ggplot(aes(x = day_1, y = day_29, color = class, shape = height)) +
  geom_point(size = 2) +
  geom_text(aes(label = gene), check_overlap = T, vjust = 1.5) +
  ylim(0, 0.03) +
  xlim(0, 0.03) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_fill_manual(values = darkpalette) +
  labs(x = "Day 1", y = "Day 29", color = "Antibiotic Class", shape = "Height") +
  theme_minimal()

