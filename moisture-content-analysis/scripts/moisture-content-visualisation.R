library(tidyverse)
library(viridis)
library(hrbrthemes)

# import moisture data from folder
moisture_data <- readxl::read_excel("moisture-content-analysis/data/raw-data/Moisture Content Analysis.xlsx",
                                    sheet = "Modified") %>%
  janitor::clean_names()

# import dataset
moisture_data_modified <- read_csv("moisture-content-analysis/data/raw-data/Moisture Content Analysis - Modified.csv", 
                                               col_types = cols(Moisture_Content_pc = col_number()))
Grouped <- Moisture_Content_Analysis_Modified %>% 
  group_by(Date, Pile_Height, Pile_Length)

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

# add weather
library(weather)
# fill in according to
# https://www.timeanddate.com/weather/uk/york/historic?month=4&year=2022
moisture_number$weather = NA
moisture_number <- mutate(moisture_number,
                       weather = case_when(
                         str_detect(date, "04-19") ~ "day-cloudy-windy",
                         str_detect(date, "04-22") ~ "day-sunny",
                         str_detect(date, "04-26") ~ "day-cloudy",
                         str_detect(date, "04-29") ~ "day-cloudy",
                         str_detect(date, "05-03") ~ "day-cloudy",
                         str_detect(date, "05-06") ~ "day-cloudy-windy",
                         str_detect(date, "05-10") ~ "day-windy",
                         str_detect(date, "05-13") ~ "day-cloudy",
                         str_detect(date, "05-17") ~ "day-light-wind"))

# change length to location A and B.
moisture_renamed <- mutate(moisture_number,
                            pile_length = case_when(
                              str_detect(pile_length, "Quarter") ~ "Y",
                              str_detect(pile_length, "Half") ~ "X"))

moisture_renamed <- mutate(moisture_renamed, pile_height = case_when(
                             str_detect(pile_height, "Bottom") ~ "bottom",
                             str_detect(pile_height, "Middle") ~ "middle",
                             str_detect(pile_height, "Top") ~ "top"))
# simple graph
moisture_renamed %>%
ggplot() +
  geom_point(aes(x = day,
                 y = moisture_content_pc,
                 colour = pile_height,
                 size = 5)) +
  geom_line(aes(x = day,
                y = moisture_content_pc,
                colour = pile_height)) +
  geom_weather(aes(x = day,
                   y = 0.85,
                   weather = weather)) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  labs(x = "day",
       y = "% moisture content",
       colour = "height") +
  scale_size(guide = 'none') +
  facet_grid(pile_length ~ .) +
  scale_color_viridis(discrete = T) +
  theme_ipsum(base_size = 10)

# second attempt on modified

# rename Length
Grouped2 <- mutate(Grouped, Pile_Length = case_when(
  str_detect(Pile_Length, "Ha") ~ "Middle",
  str_detect(Pile_Length, "Qu") ~ "Quarter"))

# remove columns
Grouped2$Boat_Weight_g = NULL
Grouped2$Start_Weight_g = NULL
Grouped2$End_Weight_g = NULL
Grouped2$Net_Start_Weight_g = NULL
Grouped2$Net_End_Weight_g = NULL

# graph
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(Grouped2, aes(y = Moisture_Content_pc, 
                     x = Date,
                     colour = Pile_Height,
                     shape = Pile_Height)) +
  geom_point(size = 3) +
  facet_wrap(~Pile_Length, dir = "v") +
  labs(y = "% Moisture Content", colour = "Pile Height", shape = "Pile Height") +
  theme_bw(base_size = 14)
