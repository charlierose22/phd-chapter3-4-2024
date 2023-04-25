library(tidyverse)

# import moisture data from folder
moisture_data <- readxl::read_excel("Moisture Content Analysis.xlsx",
                                    sheet = "Modified") %>%
  janitor::clean_names()

# change length to location A and B.
moisture_renamed <- mutate(moisture_data,
                            pile_length = case_when(
                              str_detect(sample, "Quarter") ~ "B",
                              str_detect(sample, "Half") ~ "A"))

# simple graph
moisture_renamed %>%
ggplot(aes(x = date, y = moisture_content_pc,
                         shape = pile_length,
           colour = pile_height)) +
  geom_point(size = 5) +
  geom_line() +
  theme_light(base_size = 20) +
  labs(x = "Time",
       y = "% Moisture Content",
       colour = "Height",
       shape = "Length")
