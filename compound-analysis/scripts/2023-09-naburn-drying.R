# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)
library(fuzzyjoin)
library(stringi)
library(ggforce)

# IMPORT YOUR CD DATA
drying_compounds <- readr::read_delim(
  "compound-analysis/data/raw-data/2023-09-naburn-drying.csv",
  delim = "\t",
  trim_ws = TRUE
) %>%
  janitor::clean_names()

# RENAME DATAFRAME FOR CODE TO WORK WITH MINIMAL CHANGES
basedata <- drying_compounds

# DROP THESE COLUMNS UNLESS THEY HAVE BEEN USED IN YOUR CD WORKFLOW
basedata$tags = NULL
basedata$checked = NULL

# FILTER SAMPLES WITH NO COMPOUND NAMES AND NO MS2 DATA (IF NEEDED)
basedatawithcompoundnames <- with(basedata, basedata[!(name == "" |
                                                         is.na(name)),])
ms2dataonly <- basedatawithcompoundnames[!grepl('No MS2',
                                                basedatawithcompoundnames$ms2), ]

# CREATE A UNIQUE IDENTIFIER FOR EACH FEATURE USING CONCATENATION
# CHANGE BASEDATAWITHCOMPOUNDNAMES TO MS2DATAONLY IF YOU CHOOSE TO RUN THAT CODE LINE
uniqueid <-
  add_column(basedatawithcompoundnames,
             unique_id = NA,
             .after = 0)
uniqueid$unique_id <- str_c(uniqueid$name, "_", uniqueid$rt_min)

# PEAK NUMBERS CAN BE USED AS ANOTHER IDENTIFIER
peaknumber <- add_column(uniqueid, peak_number = NA, .after = 0)
peaknumber$peak_number <- seq.int(nrow(peaknumber))

# REMOVE CD FILE NUMBERS FROM THE END OF SAMPLE NAMES
colnames(peaknumber) <-
  sub("*_raw_f\\d\\d*", "", colnames(peaknumber))

# LENGTHEN THE TABLE TO REMOVE WHITESPACE
longer <- peaknumber %>%
  pivot_longer(
    cols = group_area_a1_a:peak_rating_qc_p,
    names_to = "sample",
    values_to = "result"
  )

# CREATE A SAMPLE NAME COLUMN AND FILL, SO WE CAN GROUP PEAK RATING AND GROUP AREA
samplenames <- add_column(longer, measurement = NA)
samplenames <- mutate(samplenames,
                      measurement = case_when(
                        str_detect(sample,
                                   "group_area") ~
                          "group_area",
                        str_detect(sample,
                                   "peak_rating") ~
                          "peak_rating"
                      ))

# CLEAN SAMPLE NAME COLUMN
samplenames$sample <-
  str_replace_all(samplenames$sample, "group_area_", "")
samplenames$sample <-
  str_replace_all(samplenames$sample, "peak_rating_", "")

# WIDEN TABLE
wider <- samplenames %>%
  pivot_wider(names_from = measurement, values_from = result)

# REMOVE NAS
nona <- drop_na(wider, group_area)

# FILTER DEPENDING ON PEAK RATING NUMBER
peakrating <- subset(nona, peak_rating > 5)

# FILTER DEPENDING ON INTENSITY
grouparea <- subset(peakrating, group_area > 100000)

# FOR TECHNICAL REPLICATES ----
# CREATE A NEW COLUMN
replicates <- add_column(grouparea, replicate = NA)

# MAKE SURE TECHNICAL REPLICATES ARE AT THE END OF THE SAMPLE NAME AND CHANGE A/B/C ACCORDINGLY
replicates <- mutate(
  replicates,
  replicate = case_when(
    str_ends(sample, "a") ~ "a",
    str_ends(sample, "b") ~ "b",
    str_ends(sample, "c") ~ "c",
    str_ends(sample, "d") ~ "d",
    str_ends(sample, "e") ~ "e",
    str_ends(sample, "f") ~ "f",
    str_ends(sample, "g") ~ "g",
    str_ends(sample, "h") ~ "h",
    str_ends(sample, "i") ~ "i",
    str_ends(sample, "j") ~ "j",
    str_ends(sample, "k") ~ "k",
    str_ends(sample, "l") ~ "l",
    str_ends(sample, "m") ~ "m",
    str_ends(sample, "n") ~ "n",
    str_ends(sample, "o") ~ "o",
    str_ends(sample, "p") ~ "p",
  )
)

# ADD A SAMPLE LOCATION COLUMN AND CLEAN TO REMOVE REPLICATE NAMES SO WE CAN REMOVE SOLOS
replicates$sample_location = replicates$sample
replicates$sample_location <-
  gsub('_.*', '', replicates$sample_location)

# REMOVE PEAKS WITH RESULTS IN ONLY ONE REPLICATE
soloremoved <-
  plyr::ddply(replicates, c("unique_id", "sample_location"),
              function(d) {
                if (nrow(d) > 1)
                  d
                else
                  NULL
              })

# INCLUDE SAMPLE INFO IF NEEDED
sample_info <-
  readr::read_csv("compound-analysis/data/sample-information/2023-09-naburn-drying-samples.csv")
sample_included <- soloremoved %>% left_join(sample_info,
                                             by = "sample_location")
#sample_included$sample_name <- paste(sample_included$day,
#                                     sample_included$height,
#                                     sample_included$length,
#                                     sep = "-")
#sample_included <- select(sample_included,
#                          -day,
#                          -height,
#                          -length)

#----
# SPLIT RESULTS BASED ON MASS LIST VS MZCLOUD
split <-
  split(sample_included,
        sample_included$annot_source_mass_list_search)
mzcloud <- split$"No results"
masslists <- split$"Full match"

# MERGE MASS LISTS INTO ONE COLUMN
masslistmerged <- masslists %>%
  pivot_longer(
    cols = c(starts_with("mass_list_match")) ,
    names_to = "mass_list_name",
    names_prefix = "mass_list_match_",
    values_to = "mass_list_match"
  )

# FILTER FOR NO MATCHES AND INVALID MASS RESULTS
filteredmasslist <- masslistmerged[!grepl('No matches found',
                                          masslistmerged$mass_list_match),]
filteredmzcloud <- mzcloud[!grepl('Invalid Mass',
                                  mzcloud$annot_source_mz_cloud_search),]

# SPLIT THE INDIVIDUAL MASS LISTS
# FIRST FIND MASS LIST NAMES (APPEARING IN CONSOLE)
unique(filteredmasslist$mass_list_name)
# THEN SPLIT BY MASS LIST
splitmasslist <- split(filteredmasslist, filteredmasslist$mass_list_name)
antibiotics <- splitmasslist$"antibiotics_itn_msca_answer_160616_w_dtxsi_ds"
metabolites <- splitmasslist$"itnantibiotic_cyp_metabolites"
psychoactive <- splitmasslist$"kps_psychoactive_substances_v2"
pharmaceuticals <- splitmasslist$"kps_pharmaceuticals"

# clean compound names in antibiotics dataset.
antibiotics$name <- antibiotics$name %>%
  fedmatch::clean_strings()

# wide view for samples and fully annotated view
antibiotics_means <- antibiotics %>%
  group_by(pick(name, day, height)) %>%
  summarise(
    mean = mean(group_area),
    std = sd(group_area),
    n = length(group_area),
    se = std / sqrt(n)
  )

# add antibiotic classes for common antibiotics
class_info <-
  read.csv("compound-analysis/data/antimicrobial_classes.csv")

location_classes_means <- antibiotics_means %>%
  fuzzy_left_join(class_info,
                  by = c("name" = "name"),
                  match_fun = str_detect)

location_classes_box <- antibiotics %>%
  fuzzy_left_join(class_info,
                  by = c("name" = "name"),
                  match_fun = str_detect)

# replace NA with 'unknown'
location_classes_means$class <- location_classes_means$class %>%
  replace_na('unknown')
location_classes_box$class <- location_classes_box$class %>%
  replace_na('unknown')

# create separate data sets for location study and time series.
chem_location_study <- location_classes_means[grepl('\\<01\\>|\\<29\\>', 
                                                    location_classes_means$day), ]
chem_time_study <- location_classes_means[grepl('bottom',
                                                location_classes_means$height), ]
box_chem_location_study <- location_classes_box[grepl('\\<01\\>|\\<29\\>', 
                                                      location_classes_box$day), ]
box_chem_time_study <- location_classes_box[grepl('bottom',
                                                  location_classes_box$height), ]

chem_time_study$day = as.numeric(chem_time_study$day)
box_chem_time_study$day = as.numeric(box_chem_time_study$day)

chem_mean_time_total <- box_chem_time_study %>% 
  group_by(day, class) %>% 
  summarise(
    mean = mean(group_area),
    std = sd(group_area),
    n = length(group_area),
    se = std / sqrt(n)
  )

# PRODUCE A CSV OF RESULTS
write.csv(
  antibiotics,
  "compound-analysis/data/processed-data/2023-naburn-drying/itn_antibiotics.csv",
  row.names = FALSE
)
write.csv(
  location_classes,
  "compound-analysis/data/processed-data/2023-naburn-drying/antibiotics_class_organised.csv",
  row.names = FALSE
)
write.csv(
  metabolites,
  "compound-analysis/data/processed-data/2023-naburn-drying/itn_metabolites.csv",
  row.names = FALSE
)
write.csv(
  psychoactive,
  "compound-analysis/data/processed-data/2023-naburn-drying/psychoactive.csv",
  row.names = FALSE
)
write.csv(
  pharmaceuticals,
  "compound-analysis/data/processed-data/2023-naburn-drying/pharmaceuticals.csv",
  row.names = FALSE
)

# CHANGE HEIGHT AND WIDTH AND MIDPOINT AS NEEDED
library(RColorBrewer)
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

# Individual Plots
# split based on target antibiotics for location 
# Create a list of target antibiotics.
class_info2 <- unique(class_info$class)

# split based on target antibiotics for location 
split_location_compound <- split(chem_location_study, chem_location_study$class)
chem_loc_amino <- split_location_compound$aminoglycoside
chem_loc_fung <- split_location_compound$antifungal
chem_loc_beta <- split_location_compound$`beta-lactam`
chem_loc_glyc <- split_location_compound$glycopeptide_metronidazole
chem_loc_mac <- split_location_compound$macrolide_lincosamide
chem_loc_other <- split_location_compound$other
chem_loc_phen <- split_location_compound$phenicol
chem_loc_quin <- split_location_compound$quinolone
chem_loc_sulf <- split_location_compound$sulfonamide_trimethoprim
chem_loc_tet <- split_location_compound$tetracycline
chem_loc_unknown <- split_location_compound$unknown

split_time_compound <- split(chem_time_study, chem_time_study$class)
chem_time_amino <- split_time_compound$aminoglycoside
chem_time_fung <- split_time_compound$antifungal
chem_time_beta <- split_time_compound$`beta-lactam`
chem_time_glyc <- split_time_compound$glycopeptide_metronidazole
chem_time_mac <- split_time_compound$macrolide_lincosamide
chem_time_other <- split_time_compound$other
chem_time_phen <- split_time_compound$phenicol
chem_time_quin <- split_time_compound$quinolone
chem_time_sulf <- split_time_compound$sulfonamide_trimethoprim
chem_time_tet <- split_time_compound$tetracycline
chem_time_unknown <- split_time_compound$unknown

box_split_location_compound <- split(box_chem_location_study, box_chem_location_study$class)
box_chem_loc_amino <- box_split_location_compound$aminoglycoside
box_chem_loc_fung <- box_split_location_compound$antifungal
box_chem_loc_beta <- box_split_location_compound$`beta-lactam`
box_chem_loc_glyc <- box_split_location_compound$glycopeptide_metronidazole
box_chem_loc_mac <- box_split_location_compound$macrolide_lincosamide
box_chem_loc_other <- box_split_location_compound$other
box_chem_loc_phen <- box_split_location_compound$phenicol
box_chem_loc_quin <- box_split_location_compound$quinolone
box_chem_loc_sulf <- box_split_location_compound$sulfonamide_trimethoprim
box_chem_loc_tet <- box_split_location_compound$tetracycline
box_chem_loc_unknown <- box_split_location_compound$unknown

box_split_time_compound <- split(box_chem_time_study, box_chem_time_study$class)
box_chem_time_amino <- box_split_time_compound$aminoglycoside
box_chem_time_fung <- box_split_time_compound$antifungal
box_chem_time_beta <- box_split_time_compound$`beta-lactam`
box_chem_time_glyc <-box_split_time_compound$glycopeptide_metronidazole
box_chem_time_mac <- box_split_time_compound$macrolide_lincosamide
box_chem_time_other <- box_split_time_compound$other
box_chem_time_phen <- box_split_time_compound$phenicol
box_chem_time_quin <- box_split_time_compound$quinolone
box_chem_time_sulf <- box_split_time_compound$sulfonamide_trimethoprim
box_chem_time_tet <- box_split_time_compound$tetracycline
box_chem_time_unknown <- box_split_time_compound$unknown

# count number of genes per class
chem_count_loc <- chem_location_study %>%
  group_by(class, day, height) %>% 
  summarise(count = n_distinct(name.x))

chem_count_time <- chem_time_study %>%
  group_by(class, day) %>% 
  summarise(count = n_distinct(name.x))

chem_count_loc %>% 
  ggplot(aes(fill = class, y = count, x = day)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 10),
                    labels = c("aminoglycoside",
                               "antifungal",
                               "beta-lactam",
                               "macrolide lincosamide",
                               "other",
                               "phenicol",
                               "quinolone",
                               "sulfonamide and\ntrimethoprim",
                               "tetracycline",
                               "unknown")) +
  labs(y = "number of compounds detected", 
       fill = "antibiotic class") +
  facet_wrap(~height) +
  theme_minimal(base_size = 12)

box_chem_location_study %>%
  ggplot(aes(x = day, y = group_area, fill = class)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Spectral", n = 10),
                    labels = c("aminoglycoside",
                               "antifungal",
                               "beta-lactam",
                               "macrolide lincosamide",
                               "other",
                               "phenicol",
                               "quinolone",
                               "sulfonamide and\ntrimethoprim",
                               "tetracycline",
                               "unknown")) +
  labs(x = "day", y = "intensity", fill = "antibiotic class") +
  facet_grid(rows = vars(height)) +
  theme_minimal(base_size = 12)

# aminoglycoside location
box_chem_loc_amino %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_amino %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# beta-lac location
box_chem_loc_beta %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_beta %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# antifungal location
box_chem_loc_fung %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_fung %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# mac location
box_chem_loc_mac %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_mac %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# other location
box_chem_loc_other %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_other %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# phenicol location
box_chem_loc_phen %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_phen %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# quinolone location
box_chem_loc_quin %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_quin %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# sulf location
box_chem_loc_sulf %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_sulf %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x, ncol = 3) +
  theme_minimal(base_size = 12)

# tet location
box_chem_loc_tet %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_tet %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# unknown location
box_chem_loc_unknown %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "antiquewhite4", shape = 4) +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "height", y = "compound intensity", fill = "height") +
  facet_wrap(~day) +
  theme_minimal(base_size = 12)

# per gene
box_chem_loc_unknown %>%
  ggplot(aes(x = day, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans='log10') +
  scale_fill_manual(values = brewer.pal("Accent", n = 3)) +
  labs(x = "day", y = "relative abundance", fill = "height") +
  facet_wrap(~name.x) +
  theme_minimal(base_size = 12)

# TIME ------------
chem_mean_time_total %>%
  ggplot(aes(x = day, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "compound intensity", color = "antibiotic class") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette,
                     labels = c("aminoglycoside",
                                "antifungal",
                                "beta-lactam",
                                "macrolide lincosamide",
                                "other",
                                "phenicol",
                                "quinolone",
                                "sulfonamide and\ntrimethoprim",
                                "tetracycline",
                                "unknown")) +
  theme_minimal(base_size = 12) +
  theme(legend.position="bottom")

# aminolgycoside over time
chem_time_amino %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# beta over time
chem_time_beta %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# antifungal over time
chem_time_fung %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# mac over time
chem_time_mac %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# other over time
chem_time_other %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# phen over time
chem_time_phen %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# quin over time
chem_time_quin %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# sulf over time
chem_time_sulf %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12) +
  theme(legend.position="bottom")

# tet over time
chem_time_tet %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values = darkpalette) +
  theme_minimal(base_size = 12)

# unknown over time
chem_time_unknown %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "day", y = "intensity", color = "compound name") +
  scale_y_continuous(trans='log10') +
  theme_minimal(base_size = 12) +
  theme(legend.position="bottom")