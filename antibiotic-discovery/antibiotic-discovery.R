library(grafify)

# import dataset
discovery_timeline <- readxl::read_excel("antibiotic-discovery/antibiotic-discovery.xlsx", 
                                         col_types = c("text", "date"))

# mutate year column to years with lubridate
discovery_timeline$year <- lubridate::as_date(discovery_timeline$year)

# plot
discovery_timeline %>%
  ggplot(aes(x = year, y = 0, colour = class)) +
  geom_point(shape = 15) +
  geom_line() +
  labs(x = "Year", y = "", colour = "Antibiotic Class") +
  scale_x_date(date_breaks = "5 years", date_labels = "%y") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
