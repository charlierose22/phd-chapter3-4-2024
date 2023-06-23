library(tidyverse)
library(viridis)
library(hrbrthemes)

# import moisture data from folder
moisture_data <- readxl::read_excel("Moisture Content Analysis.xlsx",
                                    sheet = "Modified") %>%
  janitor::clean_names()

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
