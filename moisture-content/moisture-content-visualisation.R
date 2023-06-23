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
                              str_detect(pile_length, "Quarter") ~ "Z",
                              str_detect(pile_length, "Half") ~ "Y"))

# simple graph
moisture_renamed %>%
ggplot() +
  geom_point(aes(x = day,
                 y = moisture_content_pc,
                 shape = pile_length,
                 colour = pile_height,
                 size = 5)) +
  geom_line(aes(x = day,
                y = moisture_content_pc,
                shape = pile_length,
                colour = pile_height)) +
  geom_weather(aes(x = day,
                   y = 0.85,
                   weather = weather)) +
  labs(x = "day",
       y = "% moisture content",
       colour = "height",
       shape = "location") +
  scale_size(guide = 'none') +
  scale_color_viridis(discrete = T) +
  theme_ipsum(base_size = 10)
