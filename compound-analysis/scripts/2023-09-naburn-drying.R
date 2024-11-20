# LOAD TIDYVERSE AND CHECK PACKAGE UPDATES
library(tidyverse)
library(fuzzyjoin)
library(stringi)
library(ggforce)
library(RColorBrewer)
library(grafify)
library(ggrepel)
library(ggtext)
library(dygraphs)
library(emmeans)

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

# drop_na
location_classes_means$class <- location_classes_means$class %>%
  replace_na('unknown')
location_classes_box$class <- location_classes_box$class %>%
  replace_na('unknown')

location_classes_means_drop <- location_classes_means[!grepl('unknown|antifungal', 
                                                            location_classes_means$class), ]
location_classes_box_drop <- location_classes_box[!grepl('unknown|antifungal', 
                                                            location_classes_box$class), ]

# capitalise names
location_classes_means_drop$name.x <- str_to_title(location_classes_means_drop$name.x)
location_classes_box_drop$name.x <- str_to_title(location_classes_box_drop$name.x)

# create separate data sets for location study and time series.
chem_location_study <- location_classes_means_drop[grepl('\\<29\\>', 
                                                    location_classes_means_drop$day), ]
chem_time_study <- location_classes_means_drop[grepl('bottom',
                                                location_classes_means_drop$height), ]
box_chem_location_study <- location_classes_box_drop[grepl('\\<29\\>', 
                                                      location_classes_box_drop$day), ]
box_chem_time_study <- location_classes_box_drop[grepl('bottom',
                                                  location_classes_box_drop$height), ]

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
chem_loc_beta <- split_location_compound$`beta-lactam`
chem_loc_glyc <- split_location_compound$glycopeptide_metronidazole
chem_loc_mac <- split_location_compound$macrolide_lincosamide
chem_loc_other <- split_location_compound$other
chem_loc_phen <- split_location_compound$phenicol
chem_loc_quin <- split_location_compound$quinolone
chem_loc_sulf <- split_location_compound$sulfonamide_trimethoprim
chem_loc_tet <- split_location_compound$tetracycline

split_time_compound <- split(chem_time_study, chem_time_study$class)
chem_time_amino <- split_time_compound$aminoglycoside
chem_time_beta <- split_time_compound$`beta-lactam`
chem_time_glyc <- split_time_compound$glycopeptide_metronidazole
chem_time_mac <- split_time_compound$macrolide_lincosamide
chem_time_other <- split_time_compound$other
chem_time_phen <- split_time_compound$phenicol
chem_time_quin <- split_time_compound$quinolone
chem_time_sulf <- split_time_compound$sulfonamide_trimethoprim
chem_time_tet <- split_time_compound$tetracycline

box_split_location_compound <- split(box_chem_location_study, box_chem_location_study$class)
box_chem_loc_amino <- box_split_location_compound$aminoglycoside
box_chem_loc_beta <- box_split_location_compound$`beta-lactam`
box_chem_loc_glyc <- box_split_location_compound$glycopeptide_metronidazole
box_chem_loc_mac <- box_split_location_compound$macrolide_lincosamide
box_chem_loc_other <- box_split_location_compound$other
box_chem_loc_phen <- box_split_location_compound$phenicol
box_chem_loc_quin <- box_split_location_compound$quinolone
box_chem_loc_sulf <- box_split_location_compound$sulfonamide_trimethoprim
box_chem_loc_tet <- box_split_location_compound$tetracycline

box_split_time_compound <- split(box_chem_time_study, box_chem_time_study$class)
box_chem_time_amino <- box_split_time_compound$aminoglycoside
box_chem_time_beta <- box_split_time_compound$`beta-lactam`
box_chem_time_glyc <-box_split_time_compound$glycopeptide_metronidazole
box_chem_time_mac <- box_split_time_compound$macrolide_lincosamide
box_chem_time_other <- box_split_time_compound$other
box_chem_time_phen <- box_split_time_compound$phenicol
box_chem_time_quin <- box_split_time_compound$quinolone
box_chem_time_sulf <- box_split_time_compound$sulfonamide_trimethoprim
box_chem_time_tet <- box_split_time_compound$tetracycline

# edit class names for facets
box_chem_time_study$class <- str_to_title(box_chem_time_study$class)
box_chem_location_study$class <- str_to_title(box_chem_location_study$class)
chem_time_study$class <- str_to_title(chem_time_study$class)
chem_location_study$class <- str_to_title(chem_location_study$class)

box_chem_time_study$class <- str_replace_all(box_chem_time_study$class,'Beta-Lactam','Beta-lactam')
box_chem_time_study$class <- str_replace_all(box_chem_time_study$class,'Macrolide_lincosamide','Macrolides and Lincosamides')
box_chem_time_study$class <- str_replace_all(box_chem_time_study$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')
box_chem_location_study$class <- str_replace_all(box_chem_location_study$class,'Beta-Lactam','Beta-lactam')
box_chem_location_study$class <- str_replace_all(box_chem_location_study$class,'Macrolide_lincosamide','Macrolides and Lincosamides')
box_chem_location_study$class <- str_replace_all(box_chem_location_study$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')
chem_time_study$class <- str_replace_all(chem_time_study$class,'Beta-Lactam','Beta-lactam')
chem_time_study$class <- str_replace_all(chem_time_study$class,'Macrolide_lincosamide','Macrolides and Lincosamides')
chem_time_study$class <- str_replace_all(chem_time_study$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')
chem_location_study$class <- str_replace_all(chem_location_study$class,'Beta-Lactam','Beta-lactam')
chem_location_study$class <- str_replace_all(chem_location_study$class,'Macrolide_lincosamide','Macrolides and Lincosamides')
chem_location_study$class <- str_replace_all(chem_location_study$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')

# count number of genes per class
chem_count_loc <- chem_location_study %>%
  group_by(class, height) %>% 
  summarise(count = n_distinct(name.x))

chem_count_time <- chem_time_study %>%
  group_by(class, day) %>% 
  summarise(count = n_distinct(name.x))

# make factor datasets for ANOVA analysis
class_names <- c("Aminoglycosides",
                 "Beta-lactams",
                 "Glycopeptides and Metronidazole",
                 "Macrolides and Lincosamides",
                 "Other",
                 "Chloramphenicol",
                 "Quinolones",
                 "Sulfonamides and Trimethoprim",
                 "Tetracyclines")

# replace days as factor for significance tests
chem_time_factor <- box_chem_time_study
chem_time_factor$day <- as.factor(chem_time_factor$day)
chem_time_factor$name.x <- as.factor(chem_time_factor$name.x)
chem_location_factor <- box_chem_location_study
chem_location_factor$day <- as.factor(chem_location_factor$day)
chem_location_factor$height <- as.factor(chem_location_factor$height)
chem_location_factor$name.x <- as.factor(chem_location_factor$name.x)

# split based on factor dfs
# split based on target antibiotics for location 
chem_split_loc_factor <- split(chem_location_factor, chem_location_factor$class)
chem_loc_factor_amino <- chem_split_loc_factor$Aminoglycoside
chem_loc_factor_beta <- chem_split_loc_factor$'Beta-lactam'
chem_loc_factor_glyc <- chem_split_loc_factor$'Glycopeptides and Metronidazole'
chem_loc_factor_mac <- chem_split_loc_factor$'Macrolides and Lincosamides'
chem_loc_factor_other <- chem_split_loc_factor$Other
chem_loc_factor_phen <- chem_split_loc_factor$Phenicol
chem_loc_factor_quin <- chem_split_loc_factor$Quinolone
chem_loc_factor_sulf <- chem_split_loc_factor$'Sulfonamides and Trimethoprim'
chem_loc_factor_tet <- chem_split_loc_factor$Tetracycline
chem_split_time_factor <- split(chem_time_factor, chem_time_factor$class)
chem_time_factor_amino <- chem_split_time_factor$Aminoglycoside
chem_time_factor_beta <- chem_split_time_factor$'Beta-lactam'
chem_time_factor_glyc <- chem_split_time_factor$'Glycopeptides and Metronidazole'
chem_time_factor_mac <- chem_split_time_factor$'Macrolides and Lincosamides'
chem_time_factor_other <- chem_split_time_factor$Other
chem_time_factor_phen <- chem_split_time_factor$Phenicol
chem_time_factor_quin <- chem_split_time_factor$Quinolone
chem_time_factor_sulf <- chem_split_time_factor$'Sulfonamides and Trimethoprim'
chem_time_factor_tet <- chem_split_time_factor$Tetracycline

# statistics
# time study
library(FSA)
# all classes
# one way anova
one_way_all_chem <- aov(group_area ~ class, data = chem_time_factor)
summary(one_way_all_chem)
# two way anova
two_way_all_chem <- aov(group_area ~ day * class, data = chem_time_factor)
summary(two_way_all_chem)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_all_chem)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor, group_area ~ class)
# post-hoc test
chem_time_factor$class <- as.factor(chem_time_factor$class)
dunnTest(group_area ~ class,
         data = chem_time_factor,
         method = "bh")

ggplot(chem_time_factor, aes(x = class, y = group_area) ) +
  geom_boxplot() +
  scale_x_discrete(name = "Antibiotic Class") +
  scale_y_continuous(name = "UHPLC Peak Intensity") +
  theme_classic()

# Aminoglycosides
# one way anova
chem_one_way_amino <- aov(group_area ~ name.x, data = chem_time_factor_amino)
summary(chem_one_way_amino)
chem_one_way_amino_day <- aov(group_area ~ day, data = chem_time_factor_amino)
summary(chem_one_way_amino_day)
# two way anova
chem_two_way_amino <- aov(group_area ~ day * name.x, data = chem_time_factor_amino)
summary(chem_two_way_amino)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_amino)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_amino, group_area ~ name.x)
kruskal.test(data = chem_time_factor_amino, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ day,
         data = chem_time_factor_amino,
         method = "bh")

# Beta-lactams
# one way anova
chem_one_way_beta <- aov(group_area ~ name.x, data = chem_time_factor_beta)
summary(chem_one_way_beta)
chem_one_way_beta_day <- aov(group_area ~ day, data = chem_time_factor_beta)
summary(chem_one_way_beta_day)
# two way anova
chem_two_way_beta <- aov(group_area ~ day * name.x, data = chem_time_factor_beta)
summary(chem_two_way_beta)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_beta_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_beta, group_area ~ name.x)
kruskal.test(data = chem_time_factor_beta, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ name.x,
         data = chem_time_factor_beta,
         method = "bh")
dunnTest(group_area ~ day,
         data = chem_time_factor_beta,
         method = "bh")

# Macrolide Lincosamide
# one way anova
chem_one_way_mac <- aov(group_area ~ name.x, data = chem_time_factor_mac)
summary(chem_one_way_mac)
chem_one_way_mac_day <- aov(group_area ~ day, data = chem_time_factor_mac)
summary(chem_one_way_mac_day)
# two way anova
chem_two_way_mac <- aov(group_area ~ day * name.x, data = chem_time_factor_mac)
summary(chem_two_way_mac)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_mac_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_mac, group_area ~ name.x)
kruskal.test(data = chem_time_factor_mac, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ name.x,
         data = chem_time_factor_mac,
         method = "bh")
dunnTest(group_area ~ day,
         data = chem_time_factor_mac,
         method = "bh")

# Chloramphenicol
# one way anova
chem_one_way_phen_day <- aov(group_area ~ day, data = chem_time_factor_phen)
summary(chem_one_way_phen_day)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_phen_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_phen, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ day,
         data = chem_time_factor_phen,
         method = "bh")

# Quinolones
# one way anova
chem_one_way_quin <- aov(group_area ~ name.x, data = chem_time_factor_quin)
summary(chem_one_way_quin)
chem_one_way_quin_day <- aov(group_area ~ day, data = chem_time_factor_quin)
summary(chem_one_way_quin_day)
# two way anova
chem_two_way_quin <- aov(group_area ~ day * name.x, data = chem_time_factor_quin)
summary(chem_two_way_quin)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_quin_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_quin, group_area ~ name.x)
kruskal.test(data = chem_time_factor_quin, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ name.x,
         data = chem_time_factor_quin,
         method = "bh")
dunnTest(group_area ~ day,
         data = chem_time_factor_quin,
         method = "bh")

# Sulfonamides and Trimethoprim
# one way anova
chem_one_way_sulf <- aov(group_area ~ name.x, data = chem_time_factor_sulf)
summary(chem_one_way_sulf)
chem_one_way_sulf_day <- aov(group_area ~ day, data = chem_time_factor_sulf)
summary(chem_one_way_sulf_day)
# two way anova
chem_two_way_sulf <- aov(group_area ~ day * name.x, data = chem_time_factor_sulf)
summary(chem_two_way_sulf)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_sulf_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_sulf, group_area ~ name.x)
kruskal.test(data = chem_time_factor_sulf, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ name.x,
         data = chem_time_factor_sulf,
         method = "bh")
dunnTest(group_area ~ day,
         data = chem_time_factor_sulf,
         method = "bh")

# Tetracyclines
# one way anova
chem_one_way_tet <- aov(group_area ~ name.x, data = chem_time_factor_tet)
summary(chem_one_way_tet)
chem_one_way_tet_day <- aov(group_area ~ day, data = chem_time_factor_tet)
summary(chem_one_way_tet_day)
# two way anova
chem_two_way_tet <- aov(group_area ~ day * name.x, data = chem_time_factor_tet)
summary(chem_two_way_tet)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_tet_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_tet, group_area ~ name.x)
kruskal.test(data = chem_time_factor_tet, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ day,
         data = chem_time_factor_tet,
         method = "bh")

# Other
# one way anova
chem_one_way_other <- aov(group_area ~ name.x, data = chem_time_factor_other)
summary(chem_one_way_other)
chem_one_way_other_day <- aov(group_area ~ day, data = chem_time_factor_other)
summary(chem_one_way_other_day)
# two way anova
chem_two_way_other <- aov(group_area ~ day * name.x, data = chem_time_factor_other)
summary(chem_two_way_other)
# check normal distribution
par(mfrow=c(2,2))
plot(chem_one_way_other_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = chem_time_factor_other, group_area ~ name.x)
kruskal.test(data = chem_time_factor_other, group_area ~ day)
# post-hoc test
dunnTest(group_area ~ name.x,
         data = chem_time_factor_other,
         method = "bh")
dunnTest(group_area ~ day,
         data = chem_time_factor_other,
         method = "bh")

# PLOTTING
chem_count_loc %>% 
  ggplot(aes(fill = class, y = count, x = height, label = count)) + 
  geom_bar(position = "stack", stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), colour = "grey10") +
  scale_fill_manual(values = brewer.pal("Spectral", n = 8)) +
  scale_x_discrete(labels = c("~20 cm from\nthe ground", "~1m from\nthe ground", "~20 cm from\nthe surface")) +
  labs(y = "Number of Different Compounds Detected", x = "Height in the Biosolid Pile",
       fill = "Antibiotic Compound Class") +
  theme_bw(base_size = 12)

box_chem_location_study %>% 
ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  scale_y_continuous(trans = 'log10') +
  scale_fill_grafify(palette = "vibrant",
                     labels = c("~20 cm from\nthe ground", "~1m from\nthe ground", "~20 cm from\nthe surface")) +
  scale_x_discrete(labels = c("~20 cm from\nthe ground", "~1m from\nthe ground", "~20 cm from\nthe surface")) +
  labs(y = "UHPLC Compound Peak Intensity (log10)", x = "Height in the Biosolid Pile",
       fill = "Height") +
  facet_wrap(~class, ncol = 3) +
  theme_bw(base_size = 12)

# aminoglycoside total location
box_chem_loc_amino %>%
  ggplot(aes(x = height, y = group_area, fill = height)) +
  geom_boxplot() +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4) +
  scale_y_continuous(trans = "log10") +
  scale_fill_grafify(palette = "vibrant",
                     labels = c("~20 cm from\nthe ground", "~20 cm from\nthe surface")) +
  scale_x_discrete(labels = c("~20 cm from\nthe ground", "~20 cm from\nthe surface")) +
  labs(x = "Height in the Biosolid Pile", y = "UHPLC Compound Peak Intensity (log10)", fill = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
box_chem_loc_amino %>%
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

# TIME ------------

# total over time
chem_mean_time_total %>%
  ggplot(aes(x = day, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly",
                      labels = c("Aminoglycosides",
                                 "Beta-Lactams",
                                 "Macrolides and\nLincosamides",
                                 "Other",
                                 "Phenicols",
                                 "Quinolones",
                                 "Sulfonamides and\nTrimethoprim",
                                 "Tetracyclines")) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# aminolgycoside over time
chem_time_amino %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# beta over time
chem_time_beta %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# mac over time
chem_time_mac %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# other over time
chem_time_other %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# phen over time
chem_time_phen %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "vibrant") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# quin over time
chem_time_quin %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# sulf over time
chem_time_sulf %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# tet over time
chem_time_tet %>%
  ggplot(aes(x = day, y = mean, colour = name.x)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "UHPLC Peak Intensity", color = "Compound Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
