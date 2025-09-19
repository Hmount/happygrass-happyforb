#### Soil temperature at Wyoming site
####
library(tidyverse)

temp <- read.csv("data/soil_moisture/soil_moisture_wy/soil_temp_master.csv")

temp <- temp %>% drop_na()

# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
day <- c(1:4,9:12,17:20,25:28,33:36,41:44,49:52,57:60)
temp <- temp %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                            !block %in% covered ~ "cntl")) %>% 
  mutate(day = case_when(block %in% day ~ "1",
                         !block %in% day ~ "2"))


# Convert integers to factor (i.e. categories)
temp$drought <- factor(temp$drought)
temp$block <- factor(temp$block)
temp$day <- factor(temp$day)

# How does soil mositure vary as a function of drought
t.test(may23 ~ drought, data = temp)
t.test(june7 ~ drought, data = temp)
t.test(june21 ~ drought, data = temp)

# Produce boxplot in ggplot syntax
temp %>% 
  ggplot(aes(x=drought, y=may23)) +
  geom_boxplot() +
  theme_classic()
temp %>% 
  ggplot(aes(x=drought, y=june7)) +
  geom_boxplot() +
  theme_classic()
temp %>% 
  ggplot(aes(x=drought, y=june21)) +
  geom_boxplot() +
  theme_classic()

rects <- data.frame(xstart = seq(1,64,4), xend = seq(5,68,4), day = c("1","2"))

#pivot longer to look across three sampling events
temp_long <- temp %>%
  pivot_longer(
    cols = c(may23, june7, june21),   # columns to pivot
    names_to = "sampling_event",      # new column name for event
    values_to = "temp"         # new column name for temp
  )

head(temp_long)

# are temps different?
t.test(temp ~ drought, data = temp_long)
16.00833 - #drt
  15.11458 #cntl 
# = 0.89 degrees C hotter in the drought plots
