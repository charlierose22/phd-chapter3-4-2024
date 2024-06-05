library(RColorBrewer)
library(ggstatsplot)
library(ggrepel)
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

# total location
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
  labs(x = "day", y = "relative gene abundance", fill = "target antibiotic class") +
  facet_grid(rows = vars(height)) +
  theme_minimal(base_size = 12)

# location by mechanism
location_study %>%
  ggplot(aes(x = day, y = delta_ct, fill = mechanism)) +
  geom_boxplot() +
  scale_fill_manual(values = brewer.pal("Spectral", n = 9)) +
  labs(x = "day", y = "relative abundance", fill = "resistance mechanism") +
  facet_grid(rows = vars(height)) +
  theme_minimal(base_size = 12)

# ggstatsplot
location_study %>%
  grouped_ggbetweenstats(
    ## arguments relevant for ggbetweenstats
    x = day,
    y = delta_ct,
    grouping.var = class,
    ylab = "relative gene abundance",
    pairwise.display = "significant", ## display only significant pairwise comparisons
    p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
    ## arguments relevant for combine_plots
    plotgrid.args = list(nrow = 4)
  )

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
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_amino %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# beta-lactam total location
loc_box_beta %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_beta %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# glycopeptide total location
loc_box_glyc %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_glyc %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# macrolide total location
loc_box_mac %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_mac %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# mdr total location
loc_box_mdr %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_mdr %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# mge total location
loc_box_mge %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_mge %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# other total location
loc_box_other %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_other %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# phenicol total location
loc_box_phen %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_phen %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# quinolone total location
loc_box_quin %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_quin %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# sulf total location
loc_box_sulf %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_sulf %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# tet total location
loc_box_tet %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "relative abundance", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
loc_box_tet %>%
  ggplot(aes(x = day, y = delta_ct, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~gene) +
  theme_minimal(base_size = 12)

# TIME -------------------------------------
labels <- c("aminoglycoside",
           "beta-lactam",
           "glycopeptide and\nmetronidazole",
           "macrolide lincosamide",
           "multi-drug resistance",
           "mobile genetic elements\nand integrons",
           "other",
           "phenicol",
           "quinolone",
           "sulfonamide and\ntrimethoprim",
           "tetracycline")

mean_time_total %>%
  ggplot(aes(x = day, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative gene abundance") +
  scale_color_manual(values = darkpalette) +
  geom_text_repel(
    data = mean_time_total %>% filter(day == 29),
    aes(color = class, label = class),
    direction = "y",
    xlim = c(40, NA),
    hjust = 0,
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dotted",
    box.padding = .4,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 40), 
    breaks = seq(0, 30, by = 5)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

mechanism_time_total %>%
  ggplot(aes(x = day, y = mean, colour = mechanism)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative gene abundance") +
  scale_color_manual(values = darkpalette) +
  geom_text_repel(
    data = mechanism_time_total %>% filter(day == 29),
    aes(color = mechanism, label = mechanism),
    direction = "y",
    xlim = c(40, NA),
    hjust = 0,
    segment.size = .7,
    segment.alpha = .5,
    segment.linetype = "dotted",
    box.padding = .4,
    segment.curvature = -0.1,
    segment.ncp = 3,
    segment.angle = 20) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, 40), 
    breaks = seq(0, 30, by = 5)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

# aminolgycoside over time
time_amino %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# beta-lactam over time
time_beta %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# glyc over time
time_glyc %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# macrolide over time
time_mac %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# mdr over time
time_mdr %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# mge over time
time_mge %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# other over time
time_other %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# phen over time
time_phen %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# quin over time
time_quin %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# sulf over time
time_sulf %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# tet over time
time_tet %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "relative abundance", color = "gene") +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)
