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
# LOCATION ----------------------------------------------------------------

# Create a list of target antibiotics.
class <- unique(assay_samples$class)

# split based on target antibiotics for location 
split_location <- split(means_loc, means_loc$class)
loc_amino <- split_location$aminoglycoside
loc_beta <- split_location$`beta-lactam`
loc_glyc <- split_location$glycopeptide_metronidazole
loc_mac <- split_location$macrolide_lincosamide
loc_mdr <- split_location$mdr
loc_mge <- split_location$mge_integrons
loc_other <- split_location$other
loc_phen <- split_location$phenicol
loc_quin <- split_location$quinolone
loc_sulf <- split_location$sulfonamide_trimethoprim
loc_tet <- split_location$tetracycline

# split based on target antibiotics for location 
split_loc_box <- split(location_study, location_study$class)
loc_box_amino <- split_loc_box$aminoglycoside
loc_box_beta <- split_loc_box$`beta-lactam`
loc_box_glyc <- split_loc_box$glycopeptide_metronidazole
loc_box_mac <- split_loc_box$macrolide_lincosamide
loc_box_mdr <- split_loc_box$mdr
loc_box_mge <- split_loc_box$mge_integrons
loc_box_other <- split_loc_box$other
loc_box_phen <- split_loc_box$phenicol
loc_box_quin <- split_loc_box$quinolone
loc_box_sulf <- split_loc_box$sulfonamide_trimethoprim
loc_box_tet <- split_loc_box$tetracycline

# split based on target antibiotics for time 
split_time <- split(means_time, means_time$class)
time_amino <- split_time$aminoglycoside
time_beta <- split_time$`beta-lactam`
time_glyc <- split_time$glycopeptide_metronidazole
time_mac <- split_time$macrolide_lincosamide
time_mdr <- split_time$mdr
time_mge <- split_time$mge_integrons
time_other <- split_time$other
time_phen <- split_time$phenicol
time_quin <- split_time$quinolone
time_sulf <- split_time$sulfonamide_trimethoprim
time_tet <- split_time$tetracycline

split_time_box <- split(time_study, time_study$class)
time_box_amino <- split_time_box$aminoglycoside
time_box_beta <- split_time_box$`beta-lactam`
time_box_glyc <- split_time_box$glycopeptide_metronidazole
time_box_mac <- split_time_box$macrolide_lincosamide
time_box_mdr <- split_time_box$mdr
time_box_mge <- split_time_box$mge_integrons
time_box_other <- split_time_box$other
time_box_phen <- split_time_box$phenicol
time_box_quin <- split_time_box$quinolone
time_box_sulf <- split_time_box$sulfonamide_trimethoprim
time_box_tet <- split_time_box$tetracycline

# count
count_loc %>% 
  ggplot(aes(fill = class, y = count, x = day)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11),
                    labels = c("aminoglycoside",
                               "beta-lactam",
                               "glycopeptide and\nmetronidazole",
                               "macrolide lincosamide",
                               "multi-drug resistance",
                               "mobile genetic elements\nand integrons",
                               "other",
                               "phenicol",
                               "quinolone",
                               "sulfonamide and\ntrimethoprim",
                               "tetracycline")) +
  labs(y = "number of genes detected", 
       fill = "target antibiotic class") +
  facet_wrap(~height) +
  theme_minimal(base_size = 12)

location_study %>%
  ggplot(aes(x = day, y = delta_ct, fill = class)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11),
                    labels = c("aminoglycoside",
                               "beta-lactam",
                               "glycopeptide and\nmetronidazole",
                               "macrolide lincosamide",
                               "multi-drug resistance",
                               "mobile genetic elements\nand integrons",
                               "other",
                               "phenicol",
                               "quinolone",
                               "sulfonamide and\ntrimethoprim",
                               "tetracycline")) +
  labs(x = "day", y = "relative abundance", fill = "target antibiotic class") +
  facet_grid(rows = vars(height)) +
  theme_minimal(base_size = 12)

# moisture
moisture_stats %>%
  ggplot() +
  geom_point(aes(x = day,
                 y = mean,
                 colour = height)) +
  geom_line(aes(x = day,
                y = mean,
                colour = height)) +
  geom_errorbar(aes(x = day,
                    y = mean,
                    ymin = mean - sd, 
                    ymax = mean + sd,
                    colour = height), 
                width = 1) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "day",
       y = "% moisture content",
       colour = "height") +
  scale_color_manual(values = brewer.pal("Dark2", n = 4)) +
  theme_minimal(base_size = 12)

# aminoglycoside total location
loc_box_amino %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4") +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

loc_box_amino %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Pastel2", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)



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
