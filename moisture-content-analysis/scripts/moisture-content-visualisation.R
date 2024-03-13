library(tidyverse)
library(RColorBrewer)

# import moisture data from folder
moisture_data <- readxl::read_excel("moisture-content-analysis/data/raw-data/Moisture Content Analysis.xlsx",
                                    sheet = "Modified") %>%
  janitor::clean_names()

grouped <- moisture_data %>% 
  group_by(date, pile_height)

# add days instead of dates
moisture_data$day = NA
moisture_day <- mutate(moisture_data,
                        day = case_when(
                          str_detect(date, "04-19") ~ "1",
                          str_detect(date, "04-22") ~ "4",
                          str_detect(date, "04-26") ~ "8",
                          str_detect(date, "04-29") ~ "11",
                          str_detect(date, "05-03") ~ "15",
                          str_detect(date, "05-06") ~ "18",
                          str_detect(date, "05-10") ~ "22",
                          str_detect(date, "05-13") ~ "25",
                          str_detect(date, "05-17") ~ "29"))
moisture_number <- transform(moisture_day, day = as.double(day))

# rename height and length to lowercase.
moisture_renamed <- mutate(moisture_number,
                            pile_length = case_when(
                              str_detect(pile_length, "Quarter") ~ "quarter",
                              str_detect(pile_length, "Half") ~ "half"))

moisture_renamed <- mutate(moisture_renamed, pile_height = case_when(
                             str_detect(pile_height, "Bottom") ~ "bottom",
                             str_detect(pile_height, "Middle") ~ "middle",
                             str_detect(pile_height, "Top") ~ "top"))

moisture_stats <- moisture_renamed %>% 
  group_by(day, pile_height) %>% 
  summarise(mean = mean(moisture_content_pc),
            sd = sd(moisture_content_pc))

fit = aov(data = moisture_stats,
          mean ~ pile_height)

tukey <- TukeyHSD(fit)

#Tukey test representation:
plot(tukey, las = 1 , col = "forestgreen")

summary(fit)

moisture_stats %>% 
  kruskal.test(mean ~ pile_height)

# simple graph
moisture_stats %>%
ggplot() +
  geom_point(aes(x = day,
                 y = mean,
                 colour = pile_height)) +
  geom_line(aes(x = day,
                y = mean,
                colour = pile_height)) +
  geom_errorbar(aes(x = day,
                    y = mean,
                    ymin = mean - sd, 
                    ymax = mean + sd,
                    colour = pile_height), 
                width = 1) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  labs(x = "day",
       y = "% moisture content",
       colour = "height") +
  scale_color_manual(values = brewer.pal("Dark2", n = 4)) +
  theme_minimal(base_size = 12)

RColorBrewer::display.brewer.all()
