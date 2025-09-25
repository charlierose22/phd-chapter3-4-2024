library(RColorBrewer)
library(grafify)
library(ggrepel)
library(ggtext)
library(dygraphs)
library(ggtext)
library(emmeans)
library(FSA)

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

options(max.print=999999)
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
                 "Glycopeptides and Metronidazole",
                 "MLSB",
                 "Multi-Drug Resistance",
                 "MGEs and Integrons",
                 "Other",
                 "Chloramphenicol",
                 "Quinolones",
                 "Sulfonamides and Trimethoprim",
                 "Tetracyclines")

# edit class names for facets
location_study$class <- str_to_title(location_study$class)
location_study$mechanism <- str_to_title(location_study$mechanism)
means_loc$class <- str_to_title(means_loc$class)
time_study$class <- str_to_title(time_study$class)
means_time$class <- str_to_title(means_time$class)

location_study$class <- str_replace_all(location_study$class,'Beta-Lactam','Beta-lactam')
location_study$class <- str_replace_all(location_study$class,'Glycopeptide_metronidazole','Glycopeptides and Metronidazole')
location_study$class <- str_replace_all(location_study$class,'Macrolide_lincosamide','MLSB')
location_study$class <- str_replace_all(location_study$class,'Mdr','Multi-Drug Resistance')
location_study$class <- str_replace_all(location_study$class,'Mge_integrons','MGEs and Integrons')
location_study$class <- str_replace_all(location_study$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')
means_loc$class <- str_replace_all(means_loc$class,'Beta-Lactam','Beta-lactam')
means_loc$class <- str_replace_all(means_loc$class,'Glycopeptide_metronidazole','Glycopeptides and Metronidazole')
means_loc$class <- str_replace_all(means_loc$class,'Macrolide_lincosamide','MLSB')
means_loc$class <- str_replace_all(means_loc$class,'Mdr','Multi-Drug Resistance')
means_loc$class <- str_replace_all(means_loc$class,'Mge_integrons','MGEs and Integrons')
means_loc$class <- str_replace_all(means_loc$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')
time_study$class <- str_replace_all(time_study$class,'Beta-Lactam','Beta-lactam')
time_study$class <- str_replace_all(time_study$class,'Glycopeptide_metronidazole','Glycopeptides and Metronidazole')
time_study$class <- str_replace_all(time_study$class,'Macrolide_lincosamide','MLSB')
time_study$class <- str_replace_all(time_study$class,'Mdr','Multi-Drug Resistance')
time_study$class <- str_replace_all(time_study$class,'Mge_integrons','MGEs and Integrons')
time_study$class <- str_replace_all(time_study$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')
means_time$class <- str_replace_all(means_time$class,'Beta-Lactam','Beta-lactam')
means_time$class <- str_replace_all(means_time$class,'Glycopeptide_metronidazole','Glycopeptides and Metronidazole')
means_time$class <- str_replace_all(means_time$class,'Macrolide_lincosamide','MLSB')
means_time$class <- str_replace_all(means_time$class,'Mdr','Multi-Drug Resistance')
means_time$class <- str_replace_all(means_time$class,'Mge_integrons','MGEs and Integrons')
means_time$class <- str_replace_all(means_time$class,'Sulfonamide_trimethoprim','Sulfonamides and Trimethoprim')

# replace days as factor for significance tests
time_factor <- time_study
time_factor$day <- as.factor(time_factor$day)
location_factor <- location_study
location_factor$day <- as.factor(location_factor$day)
location_factor$height <- as.factor(location_factor$height)

# split based on factor dfs
# split based on target antibiotics for location 
split_loc_factor <- split(location_factor, location_factor$class)
loc_factor_amino <- split_loc_factor$Aminoglycoside
loc_factor_beta <- split_loc_factor$'Beta-lactam'
loc_factor_glyc <- split_loc_factor$'Glycopeptides and Metronidazole'
loc_factor_mac <- split_loc_factor$MLSB
loc_factor_mdr <- split_loc_factor$'Multi-Drug Resistance'
loc_factor_mge <- split_loc_factor$'MGEs and Integrons'
loc_factor_other <- split_loc_factor$Other
loc_factor_phen <- split_loc_factor$Phenicol
loc_factor_quin <- split_loc_factor$Quinolone
loc_factor_sulf <- split_loc_factor$'Sulfonamides and Trimethoprim'
loc_factor_tet <- split_loc_factor$Tetracycline
split_time_factor <- split(time_factor, time_factor$class)
time_factor_amino <- split_time_factor$Aminoglycoside
time_factor_beta <- split_time_factor$'Beta-lactam'
time_factor_glyc <- split_time_factor$'Glycopeptides and Metronidazole'
time_factor_mlsb <- split_time_factor$MLSB
time_factor_mdr <- split_time_factor$'Multi-Drug Resistance'
time_factor_mge <- split_time_factor$'MGEs and Integrons'
time_factor_other <- split_time_factor$Other
time_factor_phen <- split_time_factor$Phenicol
time_factor_quin <- split_time_factor$Quinolone
time_factor_sulf <- split_time_factor$'Sulfonamides and Trimethoprim'
time_factor_tet <- split_time_factor$Tetracycline

# statistics
# time study

# all classes
# one way anova
one_way_all_time <- aov(delta_ct ~ class, data = time_factor)
summary(one_way_all_time)
# two way anova
two_way_all_time <- aov(delta_ct ~ day * class, data = time_factor)
summary(two_way_all_time)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_all_time)
par(mfrow=c(1,1))
# anova
anova(one_way_all_time)
# tukey
emmeans(one_way_all_time, ~ class) |> pairs()
# plot
emmeans(one_way_all_time, ~ class) |> plot()
# kruskal wallis 
kruskal.test(data = time_factor, delta_ct ~ class)
kruskal.test(data = time_factor, delta_ct ~ day)
# post-hoc test
dunnTest(delta_ct ~ class,
         data = time_factor,
         method = "bh")
dunnTest(delta_ct ~ day,
         data = time_factor,
         method = "bh")

# all mechanisms
# one way anova
one_way_all_mech <- aov(delta_ct ~ mechanism, data = time_factor)
summary(one_way_all_mech)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_all_mech)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor, delta_ct ~ mechanism)
# post-hoc test
time_factor$mechanism <- as.factor(time_factor$mechanism)
dunnTest(delta_ct ~ mechanism,
         data = time_factor,
         method = "bh")

# Aminoglycosides
# one way anova
one_way_amino_gene <- aov(delta_ct ~ gene, data = time_factor_amino)
summary(one_way_amino_gene)
one_way_amino_day <- aov(delta_ct ~ day, data = time_factor_amino)
summary(one_way_amino_day)
# two way anova
two_way_amino <- aov(delta_ct ~ day * gene, data = time_factor_amino)
summary(two_way_amino)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_amino_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_amino, delta_ct ~ gene)
kruskal.test(data = time_factor_amino, delta_ct ~ day)
# post-hoc test
time_factor_amino$gene <- as.factor(time_factor_amino$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_amino,
         method = "bh")
time_factor_amino$day <- as.factor(time_factor_amino$day)
dunnTest(delta_ct ~ day,
         data = time_factor_amino,
         method = "bh")

# Beta-lactams
# one way anova
one_way_beta_gene <- aov(delta_ct ~ gene, data = time_factor_beta)
summary(one_way_beta_gene)
one_way_beta_day <- aov(delta_ct ~ day, data = time_factor_beta)
summary(one_way_beta_day)
# two way anova
two_way_beta <- aov(delta_ct ~ day * gene, data = time_factor_beta)
summary(two_way_beta)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_beta_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_beta, delta_ct ~ gene)
kruskal.test(data = time_factor_beta, delta_ct ~ day)
# post-hoc test
time_factor_beta$gene <- as.factor(time_factor_beta$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_beta,
         method = "bh")
time_factor_beta$day <- as.factor(time_factor_beta$day)
dunnTest(delta_ct ~ day,
         data = time_factor_beta,
         method = "bh")

# Glycopeptides and Metronidazole
# one way anova
one_way_glyc_gene <- aov(delta_ct ~ gene, data = time_factor_glyc)
summary(one_way_glyc_gene)
one_way_glyc_day <- aov(delta_ct ~ day, data = time_factor_glyc)
summary(one_way_glyc_day)
# two way anova
two_way_glyc <- aov(delta_ct ~ day * gene, data = time_factor_glyc)
summary(two_way_glyc)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_glyc_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_glyc, delta_ct ~ gene)
kruskal.test(data = time_factor_glyc, delta_ct ~ day)
# post-hoc test
time_factor_glyc$gene <- as.factor(time_factor_glyc$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_glyc,
         method = "bh")
time_factor_glyc$day <- as.factor(time_factor_glyc$day)
dunnTest(delta_ct ~ day,
         data = time_factor_glyc,
         method = "bh")

#MLSB
# one way anova
one_way_mlsb_gene <- aov(delta_ct ~ gene, data = time_factor_mlsb)
summary(one_way_mlsb_gene)
one_way_mlsb_day <- aov(delta_ct ~ day, data = time_factor_mlsb)
summary(one_way_mlsb_day)
# two way anova
two_way_mlsb <- aov(delta_ct ~ day * gene, data = time_factor_mlsb)
summary(two_way_mlsb)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mlsb_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_mlsb, delta_ct ~ gene)
kruskal.test(data = time_factor_mlsb, delta_ct ~ day)
# post-hoc test
time_factor_mlsb$gene <- as.factor(time_factor_mlsb$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_mlsb,
         method = "bh")
time_factor_mlsb$day <- as.factor(time_factor_mlsb$day)
dunnTest(delta_ct ~ day,
         data = time_factor_mlsb,
         method = "bh")

# Chloramphenicol
# one way anova
one_way_phen_gene <- aov(delta_ct ~ gene, data = time_factor_phen)
summary(one_way_phen_gene)
one_way_phen_day <- aov(delta_ct ~ day, data = time_factor_phen)
summary(one_way_phen_day)
# two way anova
two_way_phen <- aov(delta_ct ~ day * gene, data = time_factor_phen)
summary(two_way_phen)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_phen_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_phen, delta_ct ~ gene)
kruskal.test(data = time_factor_phen, delta_ct ~ day)
# post-hoc test
time_factor_phen$gene <- as.factor(time_factor_phen$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_phen,
         method = "bh")
time_factor_phen$day <- as.factor(time_factor_phen$day)
dunnTest(delta_ct ~ day,
         data = time_factor_phen,
         method = "bh")

# Quinolones
# one way anova
one_way_quin_gene <- aov(delta_ct ~ gene, data = time_factor_quin)
summary(one_way_quin_gene)
one_way_quin_day <- aov(delta_ct ~ day, data = time_factor_quin)
summary(one_way_quin_day)
# two way anova
two_way_quin <- aov(delta_ct ~ day * gene, data = time_factor_quin)
summary(two_way_quin)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_quin_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_quin, delta_ct ~ gene)
kruskal.test(data = time_factor_quin, delta_ct ~ day)
# post-hoc test
time_factor_quin$gene <- as.factor(time_factor_quin$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_quin,
         method = "bh")
time_factor_quin$day <- as.factor(time_factor_quin$day)
time_factor_qnrD <- time_factor_quin %>% 
  filter(gene == "qnrD")
dunnTest(delta_ct ~ day,
         data = time_factor_quin,
         method = "bh")

# Sulfonamides and Trimethoprim
# one way anova
one_way_sulf_gene <- aov(delta_ct ~ gene, data = time_factor_sulf)
summary(one_way_sulf_gene)
one_way_sulf_day <- aov(delta_ct ~ day, data = time_factor_sulf)
summary(one_way_sulf_day)
# two way anova
two_way_sulf <- aov(delta_ct ~ day * gene, data = time_factor_sulf)
summary(two_way_sulf)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_sulf_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_sulf, delta_ct ~ gene)
kruskal.test(data = time_factor_sulf, delta_ct ~ day)
# post-hoc test
time_factor_sulf$gene <- as.factor(time_factor_sulf$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_sulf,
         method = "bh")
time_factor_sulf$day <- as.factor(time_factor_sulf$day)
dunnTest(delta_ct ~ day,
         data = time_factor_sulf,
         method = "bh")

# Tetracyclines
# one way anova
one_way_tet_gene <- aov(delta_ct ~ gene, data = time_factor_tet)
summary(one_way_tet_gene)
one_way_tet_day <- aov(delta_ct ~ day, data = time_factor_tet)
summary(one_way_tet_day)
# two way anova
two_way_tet <- aov(delta_ct ~ day * gene, data = time_factor_tet)
summary(two_way_tet)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_tet_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_tet, delta_ct ~ gene)
kruskal.test(data = time_factor_tet, delta_ct ~ day)
# post-hoc test
time_factor_tet$gene <- as.factor(time_factor_tet$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_tet,
         method = "bh")
time_factor_tet$day <- as.factor(time_factor_tet$day)
dunnTest(delta_ct ~ day,
         data = time_factor_tet,
         method = "bh")

# Other
# one way anova
one_way_other_gene <- aov(delta_ct ~ gene, data = time_factor_other)
summary(one_way_other_gene)
one_way_other_day <- aov(delta_ct ~ day, data = time_factor_other)
summary(one_way_other_day)
# two way anova
two_way_other <- aov(delta_ct ~ day * gene, data = time_factor_other)
summary(two_way_other)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_other_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_other, delta_ct ~ gene)
kruskal.test(data = time_factor_other, delta_ct ~ day)
# post-hoc test
time_factor_other$gene <- as.factor(time_factor_other$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_other,
         method = "bh")
time_factor_other$day <- as.factor(time_factor_other$day)
dunnTest(delta_ct ~ day,
         data = time_factor_other,
         method = "bh")

# Multi-Drug Resistance
# one way anova
one_way_mdr_gene <- aov(delta_ct ~ gene, data = time_factor_mdr)
summary(one_way_mdr_gene)
one_way_mdr_day <- aov(delta_ct ~ day, data = time_factor_mdr)
summary(one_way_mdr_day)
# two way anova
two_way_mdr <- aov(delta_ct ~ day * gene, data = time_factor_mdr)
summary(two_way_mdr)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mdr_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_mdr, delta_ct ~ gene)
kruskal.test(data = time_factor_mdr, delta_ct ~ day)
# post-hoc test
time_factor_mdr$gene <- as.factor(time_factor_mdr$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_mdr,
         method = "bh")
time_factor_mdr$day <- as.factor(time_factor_mdr$day)
dunnTest(delta_ct ~ day,
         data = time_factor_mdr,
         method = "bh")

# MGEs and Integrons
# one way anova
one_way_mge_gene <- aov(delta_ct ~ gene, data = time_factor_mge)
summary(one_way_mge_gene)
one_way_mge_day <- aov(delta_ct ~ day, data = time_factor_mge)
summary(one_way_mge_day)
# two way anova
two_way_mge <- aov(delta_ct ~ day * gene, data = time_factor_mge)
summary(two_way_mge)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mge_day)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = time_factor_mge, delta_ct ~ gene)
kruskal.test(data = time_factor_mge, delta_ct ~ day)
# post-hoc test
time_factor_mge$gene <- as.factor(time_factor_mge$gene)
dunnTest(delta_ct ~ gene,
         data = time_factor_mge,
         method = "bh")
time_factor_mge$day <- as.factor(time_factor_mge$day)
dunnTest(delta_ct ~ day,
         data = time_factor_mge,
         method = "bh")

# location study

# all classes
# one way anova
one_way_all_loc <- aov(delta_ct ~ height, data = location_factor)
summary(one_way_all_loc)
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_all_loc)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = location_factor, delta_ct ~ class)
kruskal.test(data = location_factor, delta_ct ~ height)
kruskal.test(data = location_factor, delta_ct ~ mechanism)
# post-hoc test
location_factor$height <- as.factor(location_factor$height)
dunnTest(delta_ct ~ height,
         data = location_factor,
         method = "bh")

# aminoglycosides
# one way anova
one_way_amino <- aov(delta_ct ~ gene, data = loc_factor_amino)
summary(one_way_amino)
one_way_amino_height <- aov(delta_ct ~ height, data = loc_factor_amino)
summary(one_way_amino_height)
# two way anova
two_way_amino <- aov(delta_ct ~ height * gene, data = loc_factor_amino)
summary(two_way_amino)
TukeyHSD(two_way_amino, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_amino)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_amino, delta_ct ~ gene)
kruskal.test(data = loc_factor_amino, delta_ct ~ height)
kruskal.test(data = loc_factor_amino, delta_ct ~ day)
# post-hoc test
loc_factor_amino$height <- as.factor(loc_factor_amino$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_amino,
         method = "bh")
loc_factor_amino$gene <- as.factor(loc_factor_amino$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_amino,
         method = "bh")

# beta-lactam
# one way anova
one_way_beta <- aov(delta_ct ~ gene, data = loc_factor_beta)
summary(one_way_beta)
one_way_beta_height <- aov(delta_ct ~ height, data = loc_factor_beta)
summary(one_way_beta_height)
# two way anova
two_way_beta <- aov(delta_ct ~ height * gene, data = loc_factor_beta)
summary(two_way_beta)
TukeyHSD(two_way_beta, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_beta)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_beta, delta_ct ~ gene)
kruskal.test(data = loc_factor_beta, delta_ct ~ height)
kruskal.test(data = loc_factor_beta, delta_ct ~ day)
# post-hoc test
loc_factor_beta$height <- as.factor(loc_factor_beta$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_beta,
         method = "bh")
loc_factor_beta$gene <- as.factor(loc_factor_beta$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_beta,
         method = "bh")

# glycopeptides and metronidazole
# one way anova
one_way_glyc <- aov(delta_ct ~ gene, data = loc_factor_glyc)
summary(one_way_glyc)
one_way_glyc_height <- aov(delta_ct ~ height, data = loc_factor_glyc)
summary(one_way_glyc_height)
# two way anova
two_way_glyc <- aov(delta_ct ~ height * gene, data = loc_factor_glyc)
summary(two_way_glyc)
TukeyHSD(two_way_glyc, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_glyc)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_glyc, delta_ct ~ gene)
kruskal.test(data = loc_factor_glyc, delta_ct ~ height)
kruskal.test(data = loc_factor_glyc, delta_ct ~ day)
# post-hoc test
loc_factor_glyc$height <- as.factor(loc_factor_glyc$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_glyc,
         method = "bh")
loc_factor_glyc$gene <- as.factor(loc_factor_glyc$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_glyc,
         method = "bh")

# mlsb
# one way anova
one_way_mlsb <- aov(delta_ct ~ gene, data = loc_factor_mac)
summary(one_way_mlsb)
one_way_mlsb_height <- aov(delta_ct ~ height, data = loc_factor_mac)
summary(one_way_mlsb_height)
# two way anova
two_way_mlsb <- aov(delta_ct ~ height * gene, data = loc_factor_mac)
summary(two_way_mlsb)
TukeyHSD(two_way_mlsb, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mlsb)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_mac, delta_ct ~ gene)
kruskal.test(data = loc_factor_mac, delta_ct ~ height)
kruskal.test(data = loc_factor_mac, delta_ct ~ day)
# post-hoc test
loc_factor_mac$height <- as.factor(loc_factor_mac$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_mac,
         method = "bh")
loc_factor_mac$gene <- as.factor(loc_factor_mac$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_mac,
         method = "bh")

# phenicol
# one way anova
one_way_phen <- aov(delta_ct ~ gene, data = loc_factor_phen)
summary(one_way_phen)
one_way_phen_height <- aov(delta_ct ~ height, data = loc_factor_phen)
summary(one_way_phen_height)
# two way anova
two_way_phen <- aov(delta_ct ~ height * gene, data = loc_factor_phen)
summary(two_way_phen)
TukeyHSD(two_way_phen, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_phen)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_phen, delta_ct ~ gene)
kruskal.test(data = loc_factor_phen, delta_ct ~ height)
kruskal.test(data = loc_factor_phen, delta_ct ~ day)
# post-hoc test
loc_factor_phen$height <- as.factor(loc_factor_phen$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_phen,
         method = "bh")
loc_factor_phen$gene <- as.factor(loc_factor_phen$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_phen,
         method = "bh")

# quinolone
# one way anova
one_way_quin <- aov(delta_ct ~ gene, data = loc_factor_quin)
summary(one_way_quin)
one_way_quin_height <- aov(delta_ct ~ height, data = loc_factor_quin)
summary(one_way_quin_height)
# two way anova
two_way_quin <- aov(delta_ct ~ height * gene, data = loc_factor_quin)
summary(two_way_quin)
TukeyHSD(two_way_quin, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_quin)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_quin, delta_ct ~ gene)
kruskal.test(data = loc_factor_quin, delta_ct ~ height)
kruskal.test(data = loc_factor_quin, delta_ct ~ day)
# post-hoc test
loc_factor_quin$height <- as.factor(loc_factor_quin$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_quin,
         method = "bh")
loc_factor_quin$gene <- as.factor(loc_factor_quin$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_quin,
         method = "bh")

# sulfonamides and trimethoprim
# one way anova
one_way_sulf <- aov(delta_ct ~ gene, data = loc_factor_sulf)
summary(one_way_sulf)
one_way_sulf_height <- aov(delta_ct ~ height, data = loc_factor_sulf)
summary(one_way_sulf_height)
# two way anova
two_way_sulf <- aov(delta_ct ~ height * gene, data = loc_factor_sulf)
summary(two_way_sulf)
TukeyHSD(two_way_sulf, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_sulf)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_sulf, delta_ct ~ gene)
kruskal.test(data = loc_factor_sulf, delta_ct ~ height)
kruskal.test(data = loc_factor_sulf, delta_ct ~ day)
# post-hoc test
loc_factor_sulf$height <- as.factor(loc_factor_sulf$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_sulf,
         method = "bh")
loc_factor_sulf$gene <- as.factor(loc_factor_sulf$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_sulf,
         method = "bh")

# tetracycline
# one way anova
one_way_tet <- aov(delta_ct ~ gene, data = loc_factor_tet)
summary(one_way_tet)
one_way_tet_height <- aov(delta_ct ~ height, data = loc_factor_tet)
summary(one_way_tet_height)
# two way anova
two_way_tet <- aov(delta_ct ~ height * gene, data = loc_factor_tet)
summary(two_way_tet)
TukeyHSD(two_way_tet, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_tet)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_tet, delta_ct ~ gene)
kruskal.test(data = loc_factor_tet, delta_ct ~ height)
kruskal.test(data = loc_factor_tet, delta_ct ~ day)
# post-hoc test
loc_factor_tet$height <- as.factor(loc_factor_tet$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_tet,
         method = "bh")
loc_factor_tet$gene <- as.factor(loc_factor_tet$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_tet,
         method = "bh")

# other
# one way anova
one_way_other <- aov(delta_ct ~ gene, data = loc_factor_other)
summary(one_way_other)
one_way_other_height <- aov(delta_ct ~ height, data = loc_factor_other)
summary(one_way_other_height)
# two way anova
two_way_other <- aov(delta_ct ~ height * gene, data = loc_factor_other)
summary(two_way_other)
TukeyHSD(two_way_other, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_other)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_other, delta_ct ~ gene)
kruskal.test(data = loc_factor_other, delta_ct ~ height)
kruskal.test(data = loc_factor_other, delta_ct ~ day)
# post-hoc test
loc_factor_other$height <- as.factor(loc_factor_other$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_other,
         method = "bh")
loc_factor_other$gene <- as.factor(loc_factor_other$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_other,
         method = "bh")

# mdr
# one way anova
one_way_mdr <- aov(delta_ct ~ gene, data = loc_factor_mdr)
summary(one_way_mdr)
one_way_mdr_height <- aov(delta_ct ~ height, data = loc_factor_mdr)
summary(one_way_mdr_height)
# two way anova
two_way_mdr <- aov(delta_ct ~ height * gene, data = loc_factor_mdr)
summary(two_way_mdr)
TukeyHSD(two_way_mdr, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mdr)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_mdr, delta_ct ~ gene)
kruskal.test(data = loc_factor_mdr, delta_ct ~ height)
kruskal.test(data = loc_factor_mdr, delta_ct ~ day)
# post-hoc test
loc_factor_mdr$height <- as.factor(loc_factor_mdr$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_mdr,
         method = "bh")
loc_factor_mdr$gene <- as.factor(loc_factor_mdr$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_mdr,
         method = "bh")

# mge
# one way anova
one_way_mge <- aov(delta_ct ~ gene, data = loc_factor_mge)
summary(one_way_mge)
one_way_mge_height <- aov(delta_ct ~ height, data = loc_factor_mge)
summary(one_way_mge_height)
# two way anova
two_way_mge <- aov(delta_ct ~ height * gene, data = loc_factor_mge)
summary(two_way_mge)
TukeyHSD(two_way_mge, which = 'height:gene')
# check normal distribution
par(mfrow=c(2,2))
plot(one_way_mge)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = loc_factor_mge, delta_ct ~ gene)
kruskal.test(data = loc_factor_mge, delta_ct ~ height)
kruskal.test(data = loc_factor_mge, delta_ct ~ day)
# post-hoc test
loc_factor_mge$height <- as.factor(loc_factor_mge$height)
dunnTest(delta_ct ~ height,
         data = loc_factor_mge,
         method = "bh")
loc_factor_mge$gene <- as.factor(loc_factor_mge$gene)
dunnTest(delta_ct ~ gene,
         data = loc_factor_mge,
         method = "bh")

# plots
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
  geom_boxplot(show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~class, ncol = 3) +
  theme_bw(base_size = 12)

# location by mechanism
location_study %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~mechanism, ncol = 3) +
  theme_bw(base_size = 12)

# aminoglycoside total location
loc_box_amino %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_amino %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# beta-lactam total location
loc_box_beta %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_beta %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# glycopeptide total location
loc_box_glyc %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_glyc %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# macrolide total location
loc_box_mac %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_mac %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# mdr total location
loc_box_mdr %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_mdr %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# mge total location
loc_box_mge %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_mge %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# other total location
loc_box_other %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_other %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# phenicol total location
loc_box_phen %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_phen %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# quinolone total location
loc_box_quin %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_quin %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# sulf total location
loc_box_sulf %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_sulf %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  facet_wrap(~gene) +
  theme_bw(base_size = 12)

# tet total location
loc_box_tet %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
  theme_bw(base_size = 12)

# per gene
loc_box_tet %>%
  ggplot(aes(x = height, y = delta_ct, fill = height)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(position=position_jitter(0.2), color = "grey35", shape = 4, show.legend = FALSE) +
  scale_fill_grafify(palette = "vibrant") +
  scale_x_discrete(labels = c("~20 cm from\n ground", "~1m from\n ground", "~20 cm from\n surface")) +
  labs(y = "Relative Gene Abundance", x = "Height in the Biosolid Pile") +
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
            "Sulfonamides and\nTrimethoprim",
            "Tetracyclines")

# total gene over time
mean_time_total %>%
  ggplot(aes(x = day, y = mean, colour = class)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Average Relative Gene Abundance", color = "Antibiotic Class") +
  scale_color_grafify(palette = "kelly",
                      labels = c("Aminoglycosides",
                                 "Beta-lactams",
                                 "Glycopeptides and\nMetronidazole",
                                 "MLSB",
                                 "Multi-drug Resistance",
                                 "MGEs and Integrons",
                                 "Other",
                                 "Chloramphenicol",
                                 "Quinolones",
                                 "Sulfonamides and\nTrimethoprim",
                                 "Tetracyclines")) +
#  scale_y_continuous(trans='log10') +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# stacked area plot
stacked_time_gene <- time_study %>% 
  group_by(day, class) %>%
  summarise(n = sum(delta_ct)) %>%
  mutate(percentage = n / sum(n))

# time_gene
ggplot(stacked_time_gene, aes(x=day, y=percentage, fill=class)) + 
  geom_area(alpha=0.6 , linewidth=0.5, colour="black") +
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

# total mechanism over time
mechanism_time_total %>%
  ggplot(aes(x = day, y = mean, colour = mechanism)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Average Relative Gene Abundance", color = "Type of Resistance\nMechanism") +
  scale_color_grafify(palette = "kelly",
                      labels = c("Antibiotic Deactivation",
                                 "Antibiotic Target Alteration",
                                 "Antibiotic Target Protection",
                                 "Antibiotic Target Replacement",
                                 "Efflux Pumps",
                                 "Insertion Sequence",
                                 "Integrase",
                                 "Phage",
                                 "Transposase")) +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# aminolgycoside over time
time_amino %>%
  ggplot(aes(x = day, y = mean, colour = gene)) +
  geom_point(shape = 15) +
  geom_errorbar(aes(x = day,
                    ymin = mean - se,
                    ymax = mean + se),
                width = .6) +
  geom_line() +
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
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
  labs(x = "Day of Drying Period", y = "Relative Gene Abundance", color = "Gene Name") +
  scale_color_grafify(palette = "kelly") +
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
