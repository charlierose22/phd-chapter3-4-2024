# For manipulating data
library(tidyverse)
library(magrittr)
library(lubridate)
library(Hmisc)

# For descriptive statistics and graphing
library(skimr)
library(scales)
library(gridExtra)
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

moisture_stats_height <- moisture_number %>% 
  group_by(day, pile_height) %>% 
  summarise(mean = mean(moisture_content_pc),
            sd = sd(moisture_content_pc))
moisture_stats_height$mean <- moisture_stats_height$mean * 100

moisture_stats <- moisture_number %>% 
  group_by(day) %>% 
  summarise(mean = mean(moisture_content_pc),
            sd = sd(moisture_content_pc))
moisture_stats$mean <- moisture_stats$mean * 100

moisture_stats %>% skim()
moisture_stats_height %>% skim()

# simple graph
moisture_stats %>%
ggplot() +
  geom_point(aes(x = day,
                 y = mean)) +
  geom_line(aes(x = day,
                y = mean)) +
  geom_errorbar(aes(x = day,
                    y = mean,
                    ymin = mean - sd, 
                    ymax = mean + sd), 
                width = 1) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0,100)) +
  labs(x = "Day",
       y = "% Moisture Content") +
  theme_minimal(base_size = 12)

# simple graph
moisture_stats_height %>%
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
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(0,100)) +
  labs(x = "Day",
       y = "% Moisture Content",
       colour = "Height") +
  theme_minimal(base_size = 12)
