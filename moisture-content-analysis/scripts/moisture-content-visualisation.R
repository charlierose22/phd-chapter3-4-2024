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
library(grafify)
library(FSA)
library(rstatix)

# import moisture data from folder
moisture_data <- readxl::read_excel("moisture-content-analysis/data/raw-data/Moisture Content Analysis.xlsx",
                                    sheet = "Modified") %>%
  janitor::clean_names()

weather_data <- readxl::read_excel("moisture-content-analysis/data/Weather Data Topcliffe.xlsx") %>%
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
moisture_stats_height %>% skim

weather_data$day = row_number(weather_data)

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
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  labs(x = "Day of Drying Period",
       y = "% Moisture Content of the Biosolid Pile") +
  theme_bw(base_size = 12)

# weather graph
weather_data %>%
  ggplot(aes(x = day, y = avg_temp, ymin = min_temp, ymax = max_temp)) + 
  geom_line(color = "black", size = 1) +
  geom_ribbon(fill = "indianred", alpha = 0.3) +
  geom_bar(aes(y = precipitation), stat = "identity", size = .1, fill = "lightseagreen", alpha = .4) + 
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(name = "Temperature (°C)", 
                     sec.axis = sec_axis(~ . * 1, name = "Precipitation (mm)")) +
  labs(x = "Day of Drying Period") +
  theme_bw(base_size = 12)

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
  scale_x_continuous(limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20, 25, 30)) +
  scale_colour_grafify(palette = "vibrant",
                     labels = c("~20 cm from\nthe ground", "~1m from\nthe ground", "~20 cm from\nthe surface")) +
  labs(x = "Day of Drying Period",
       y = "% Moisture Content of the Biosolid Pile",
       colour = "Height in the\nBiosolid Pile") +
  theme_bw(base_size = 12)

# statistics time
moisture_stats <- moisture_number %>% 
  group_by(day) %>% 
  summarise(mean = mean(moisture_content_pc),
            sd = sd(moisture_content_pc))
moisture_number$moisture_content_pc <- moisture_number$moisture_content_pc * 100
moisture_stats$mean <- moisture_stats$mean * 100
summary(moisture_number)
# one way anova
one.way <- aov(moisture_content_pc ~ day, data = moisture_number)
summary(one.way)
# check normal distribution
par(mfrow=c(2,2))
plot(one.way)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = moisture_number, moisture_content_pc ~ day)
# post-hoc test
moisture_factor <- moisture_number
moisture_factor$day <- as.factor(moisture_factor$day)
dunnTest(moisture_content_pc ~ day,
         data = moisture_factor,
         method = "none")

# statistics height
moisture_stats_height <- moisture_number %>% 
  group_by(day, pile_height) %>% 
  summarise(mean = mean(moisture_content_pc),
            sd = sd(moisture_content_pc))
moisture_factor$pile_height <- as.factor(moisture_factor$pile_height)
moisture_stats_height$mean <- moisture_stats_height$mean * 100
summary(moisture_number)
# one way anova
one.way2 <- aov(moisture_content_pc ~ pile_height, data = moisture_number)
summary(one.way2)
# check normal distribution
par(mfrow=c(2,2))
plot(one.way2)
par(mfrow=c(1,1))
# kruskal wallis 
kruskal.test(data = moisture_number, moisture_content_pc ~ pile_height)
# post-hoc test
dunnTest(moisture_content_pc ~ pile_height,
         data = moisture_factor,
         method = "bh")