#### Soil Moisture figures (continuous measure probes)

library(tidyverse)

## read in data
## CA on other .r script from Jen (make_soil_moistire_figure.R)
#sm.ca <- read.csv("data/soil_moisture_ca/three_year_soil_moisture.csv", col.names = T)

## WY (all values in m3/m3 VWC)
sm.wy21 <- read.csv("data/soil_moisture_wy/Restoration_soilmoisture_2021.csv")
sm.wy21 <- sm.wy21 %>% mutate(across(-1, as.numeric))
sm.wy22 <- read.csv("data/soil_moisture_wy/Restoration_soilmoisture_2022.csv")
sm.wy22 <- sm.wy22 %>% mutate(across(-1, as.numeric))
sm.wy23 <- read.csv("data/soil_moisture_wy/Restoration_soilmoisture_2023(toAug).csv")
sm.wy23 <- sm.wy23 %>% mutate(across(-1, as.numeric))
#sm.wy24 <- read.csv("data/soil_moisture_wy/Restoration_soilmoisture_2024.csv")
sm.wy <- bind_rows(sm.wy21,sm.wy22,sm.wy23)
sm.wy <- sm.wy %>% mutate(Timestamp = mdy_hm(Measurement.Time))

#pivot longer and combine timestamps into daily sm
sm.wy.long <- sm.wy %>% pivot_longer(where(is.numeric), names_to = "probenumber", values_to = "Moisture")
sm.wy.long$Block <- gsub("X(\\d+)_Port\\d+", "\\1", sm.wy.long$probenumber)  # Extract block number
sm.wy.long$Port <- gsub(".*_Port([0-9]+)", "\\1", sm.wy.long$probenumber)   # Extract port number
# drought blocks
droughtblocks <- data.frame(Block= c(20,22,27,30,33,42,49,50,51,58),
                            trt = c("drt","drt","cntl","drt","cntl","cntl","cntl","drt","cntl","drt"))
sm.wy.long <- merge(sm.wy.long, droughtblocks)

test <- (sm.wy.long %>% mutate(Date = floor_date(Timestamp, unit = "day"))
         %>% select(Date, trt, Moisture,Block)
         %>% group_by(trt, Date)
         %>% summarize(Moisture = mean(Moisture, na.rm = TRUE))
         %>% arrange(trt, Date)
  )

smwy.fig <- ggplot(test, aes(
  x = Date,
  y = Moisture,
  color = trt
)) +
  geom_line(size = 1, alpha = 0.9) +
  ylab(expression(
    Soil ~ volumetric ~ water ~ content ~ (m ^ 3 ~ m ^ -3)
  )) +
  xlab(NULL) +
  #xlim(min(df$Date), max(df$Date) + years(1)) +
  #scale_x_continuous(breaks = pretty(df$Date, n = 5)) +
  
  scale_color_brewer(type = 'seq', palette = 'Blues') +
  theme_light(base_size = 14) +
  theme(plot.margin = unit(c(
    0.75, 0.75, 0.75, 0.75
  ), "in"))

library(ggpubr)
ggarrange(gg,smwy.fig, nrow=2)


t.test(test$Moisture~test$trt)


#### together
test$site <- "WY"
smwy <- test %>% select(trt,Date,Moisture,site)
names(smwy)[names(smwy) == 'trt'] <- 'Treatment'
smwy$Treatment <- as.factor(smwy$Treatment)
df$Treatment <- factor( df$Treatment , ordered = FALSE )
df$site <- "CA"
allsm <- bind_rows(smwy,df)

ggplot(allsm, (aes(y=Moisture, x=Date, col=Treatment)))+
  geom_line()+
  facet_wrap(~site)+
  scale_color_manual(values = c("blue","skyblue","skyblue","blue"))


#seperate and combine
# Create a data.frame for the shaded regions
wyshade_data <- data.frame(
  xmin = as.POSIXct(c("2021-05-01", "2022-05-01", "2023-05-01"), tz = "UTC"),
  xmax = as.POSIXct(c("2021-09-01", "2022-09-01", "2023-09-01"), tz = "UTC"),
  ymin = -Inf,  # Extend shading to the bottom of the plot
  ymax = Inf    # Extend shading to the top of the plot
)
wysmplot <- ggplot(smwy, (aes(y=Moisture, x=Date, col=Treatment)))+
  geom_line()+
  facet_wrap(~site)+
  # ylab(expression(
  #   Soil ~ volumetric ~ water ~ content ~ (m ^ 3 ~ m ^ -3)
  # ))+
  scale_x_datetime(
    limits = as.POSIXct(c("2021-01-01 00:00:00", "2023-09-01 00:00:00"), tz = "UTC"),
    date_labels = "%Y-%b",  # Optional: Custom labels like "2021-Jan"
    breaks = x_breaks
  )+
  scale_y_continuous(breaks=c(.10,.15,.20,.25,.30,.35))+
  labs(col="Precipitation 
Treatment", x=" ", y=" ")+
  scale_color_manual(values = c("blue","skyblue"), labels=c("Ambient", "Reduction"))+
  geom_rect(
    data = wyshade_data, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
    inherit.aes = F,  # Prevents inheriting aesthetics from the main plot
    fill = "grey18",        # Choose a shading color
    alpha = 0.2           # Set transparency
  ) +
  theme_minimal()+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),    # Reduce the size of legend keys
        legend.text = element_text(size = 8), # Make legend text smaller
        legend.title = element_text(size = 9),# Adjust title size (optional)
        legend.spacing.y = unit(0.2, "cm"))

cashade_data <- data.frame(
  xmin = as.POSIXct(c("2021-01-01", "2022-01-01", "2023-01-01"), tz = "UTC"),
  xmax = as.POSIXct(c("2021-04-01", "2022-04-01", "2023-04-01"), tz = "UTC"),
  ymin = -Inf,  # Extend shading to the bottom of the plot
  ymax = Inf    # Extend shading to the top of the plot
)
x_breaks <- seq(as.POSIXct("2021-01-01"), as.POSIXct("2023-09-01"), by = "6 months")  # Example: every 6 months
df$Treatment <- ordered(df$Treatment, levels = c("125%","50%"))
casmplot <- ggplot(df, (aes(y=Moisture, x=Date, col=Treatment)))+
  geom_line()+
  facet_wrap(~site)+
  # ylab(expression(
  #   Soil ~ volumetric ~ water ~ content ~ (m ^ 3 ~ m ^ -3)
  # ))+
  scale_x_datetime(
    limits = as.POSIXct(c("2021-01-01 00:00:00", "2023-09-01 00:00:00"), tz = "UTC"),
    date_labels = "%Y-%b",  # Optional: Custom labels like "2021-Jan"
    breaks = x_breaks
  )+
  scale_y_continuous(breaks=c(.10,.15,.20,.25,.30))+
  labs(col="Precipitation
Treatment", x=" ",y=" ")+
  scale_color_manual(values = c("blue","skyblue"), labels=c("Addition", "Reduction"))+
  geom_rect(
    data = cashade_data, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
    inherit.aes = F,  # Prevents inheriting aesthetics from the main plot
    fill = "grey18",        # Choose a shading color
    alpha = 0.2           # Set transparency
  ) +
  theme_minimal()+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),    # Reduce the size of legend keys
        legend.text = element_text(size = 8), # Make legend text smaller
        legend.title = element_text(size = 9),# Adjust title size (optional)
        legend.spacing.y = unit(0.2, "cm"))

library(ggpubr)
smwy$Treatment <- relevel(smwy$Treatment, ref = "cntl")
df$Treatment <- relevel(df$Treatment, ref = "125%")

smcombo<- ggarrange(wysmplot,casmplot, nrow=2)
smcombo<-annotate_figure(smcombo, left=text_grob(expression(
  Soil ~ volumetric ~ water ~ content ~ (m ^ 3 ~ m ^ -3)
),rot=90))
smcombo

##analysis
t.test(smwy$Moisture~smwy$Treatment)
((.1838059-.2019800)/.2018900)*100 # percent reduction CA
t.test(df$Moisture~df$Treatment)
((.1596717-.1819071)/.1819071)*100 # percent reduction CA
