library(tidyverse)
library(grafify)

# import dataset
discovery_timeline <- readxl::read_excel("antibiotic-discovery/antibiotic-discovery.xlsx", 
                                         col_types = c("text", "date", "text"))

# mutate year column to years with lubridate
discovery_timeline$year <- lubridate::as_date(discovery_timeline$year)
discovery_timeline$year2 <- lubridate::year(discovery_timeline$year)

# plot
discovery_timeline %>%
  ggplot(aes(x = year, y = 0.5, colour = class, label = class)) +
  geom_point(shape = 16, size = 4) +
  labs(x = "Year", y = "") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_discrete(labels = "") +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

