# load libraries
library(tidyverse)
library(readxl)

# import data
top_30_genus_days <- read_excel("microbial-community-analysis/novogene-edited/top-30-genus-days.xlsx")
top_30_genus_height <- read_excel("microbial-community-analysis/novogene-edited/top-30-genus-height.xlsx")

# sample annotation:
# A4 - Bottom Day 1
# A5 - Middle Day 1
# A6 - Top Day 1
# A52 - Bottom Day 29
# A53 - Middle Day 29
# A54 - Top Day 29

# pivot data
days_long <- pivot_longer(top_30_genus_days, cols = c(Methanosaeta:Others),
                          names_to = "genus", values_to = "abundance")
height_long <- pivot_longer(top_30_genus_height, cols = c(Methanosaeta:Others),
                            names_to = "genus", values_to = "abundance")

# recode day column
days_long$day = NA
days <- mutate(days_long, day = case_when(str_detect(Day, "Day 1") ~ "1",
                                          str_detect(Day, "Day 29") ~ "29"))
colnames(height_long) <- c("height", "genus", "abundance")
heights <- height_long

# plot
ggplot(days, aes(x = day, y = abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  scale_fill_grafify(palette = "kelly") +
  labs(x = "Day of Drying Period", y = "Proportion of Relative Gene Abundance", fill = "Genus of Microbes") +
  theme_bw(base_size = 12)

ggplot(heights, aes(x = height, y = abundance, fill = genus)) +
  geom_bar(stat = "identity") +
  scale_fill_grafify(palette = "kelly") +
  labs(x = "Height in Biosolid Pile", y = "Proportion of Relative Gene Abundance", fill = "Genus of Microbes") +
  theme_bw(base_size = 12)
