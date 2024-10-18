library(RColorBrewer)
library(grafify)
library(ggrepel)
library(ggtext)
library(dygraphs)
library(ggtext)
darkpalette <- c("firebrick4",
                 "orangered1",
                 "darkgoldenrod2",
                 "yellow3",
                 "darkolivegreen4", 
                 "green3",
                 "turquoise3", 
                 "dodgerblue3", 
                 "darkslateblue", 
                 "darkorchid2",
                 "violetred2",
                 "maroon")
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

class_names <- c("Aminoglycosides",
                 "Beta-lactams",
                 "Glycopeptides and metronidazole",
                 "MLSB",
                 "Multi-drug resistance",
                 "MGEs and integrons",
                 "Other",
                 "Chloramphenicol",
                 "Quinolones",
                 "Sulfonamides and trimethoprim",
                 "Tetracyclines")

# count
count_loc %>% 
  ggplot(aes(fill = class, y = count, x = height, label = count)) + 
  geom_bar(position = "stack", stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), colour = "grey10") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11),
                    labels = c("Aminoglycosides",
                               "Beta-lactams",
                               "Glycopeptides and\nMetronidazole",
                               "MLSB",
                               "Multi-drug resistance",
                               "MGEs and Integrons",
                               "Other",
                               "Phenicols",
                               "Quinolones",
                               "Sulfonamides and\nTrimethoprim",
                               "Tetracyclines")) +
  scale_x_discrete(labels = c("~20 cm from\nthe ground", "~1m from\nthe ground", "~20 cm from\nthe surface")) +
  labs(y = "Number of Different Genes Detected", x = "Height in the Biosolid Pile",
       fill = "Target Antibiotic Class") +
  theme_bw(base_size = 12)

# total location
location_study %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans = 'log10') +
  scale_fill_grafify(palette = "vibrant",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative gene abundance", x = "Pile height",
       fill = "Pile height") +
  facet_wrap(~class, ncol = 3) +
  theme_bw(base_size = 12)

# location by mechanism
location_study %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans = 'log10') +
  scale_fill_grafify(palette = "vibrant",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative gene abundance", x = "Pile height",
       fill = "Pile height") +
  facet_wrap(~mechanism, ncol = 3) +
  theme_bw(base_size = 12)

# aminoglycoside total location
loc_box_amino %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_amino %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# beta-lactam total location
loc_box_beta %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_beta %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# glycopeptide total location
loc_box_glyc %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_glyc %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# macrolide total location
loc_box_mac %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_mac %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# mdr total location
loc_box_mdr %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_mdr %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# mge total location
loc_box_mge %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_mge %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# other total location
loc_box_other %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_other %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# phenicol total location
loc_box_phen %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_phen %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# quinolone total location
loc_box_quin %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_quin %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# sulf total location
loc_box_sulf %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_sulf %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# tet total location
loc_box_tet %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  theme_bw(base_size = 12)

# per gene
loc_box_tet %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_grafify(palette = "okabe_ito",
                     labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(x = "Pile height", y = "Relative gene abundance", fill = "Pile height") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# TIME -------------------------------------
labels <- c("Aminoglycosides",
            "Beta-lactams",
            "Glycopeptides and\nMetronidazole",
            "MLSB",
            "Multi-drug Resistance",
            "MGEs and Integrons",
            "Other",
            "Chloramphenicol",
            "Quinolones",
            "Sulfonamides and\ntrimethoprim",
            "Tetracyclines")

# stacked area plot
stacked_time_gene <- time_study %>% 
  group_by(day, class) %>%
  summarise(n = sum(delta_ct)) %>%
  mutate(percentage = n / sum(n))

# time_gene
ggplot(stacked_time_gene, aes(x=day, y=percentage, fill=class)) + 
  geom_area(alpha=0.6 , size=0.5, colour="black") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 11),
                    labels = c("Aminoglycosides",
                               "Beta-Lactams",
                               "Glycopeptides and\nMetronidazole",
                               "MLSB",
                               "Multi-Drug Resistance",
                               "MGEs and Integrons",
                               "Other",
                               "Chloramphenicol",
                               "Quinolones",
                               "Sulfonamides and\nTrimethoprim",
                               "Tetracyclines")) +
  labs(x = "Day of Drying Period", y = "Proportion of Relative Gene Abundance", fill = "Target Antibiotic Class") +
  theme_bw(base_size = 12)

# stacked area plot
stacked_time_mechanism <- time_study %>% 
  group_by(day, mechanism) %>% 
  summarise(n = sum(delta_ct)) %>%
  mutate(percentage = n / sum(n))

# time_mechanism
ggplot(stacked_time_mechanism, aes(x = day, y = percentage, fill = mechanism)) + 
  geom_area(alpha = 0.6 , size = 0.5, colour = "black") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 9),
                    labels = c("Antibiotic Deactivation",
                               "Antibiotic Target Alteration",
                               "Antibiotic Target Protection",
                               "Antibiotic Target Replacement",
                               "Efflux Pumps",
                               "Insertion Sequence",
                               "Integrase",
                               "Phage",
                               "Transposase")) +
  labs(x = "Day of Drying Period", y = "Proportion of Relative Gene Abundance", fill = "Type of Resistance\nMechanism") +
  theme_bw(base_size = 12)

# aminolgycoside over time
time_amino %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# beta-lactam over time
time_beta %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# glyc over time
time_glyc %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# macrolide over time
time_mac %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# mdr over time
time_mdr %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# mge over time
time_mge %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# other over time
time_other %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# phen over time
time_phen %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# quin over time
time_quin %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# sulf over time
time_sulf %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# tet over time
time_tet %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance (∆Ct)", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
