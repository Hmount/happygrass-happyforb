#### Cleaning data, making figures, and analyzing differences in Soil Moisture 
#### measured from continuous probes.
#### First, both datasets are cleaned and prepared. Next, a nice joint figure 
#### with both sites is produced. Finally, t-tests are used to analyze 
#### significant differences by precipitation treatment and percent reduction. 
#### (note: differences in treatments during the growing season are reported).

library(tidyverse)
library(ggpubr)

## read in + clean data by site
## produce nice single-site figure for each
## CA (all values in m3/m3 VWC)
header <- (scan("data/soil_moisture/soil_moisture_ca/three_year_soil_moisture.csv",
                nlines = 1,
                sep=",",
                what = character())
           %>% tail(-1))

sm.ca <- (
  readr::read_csv('data/soil_moisture/soil_moisture_ca/three_year_soil_moisture.csv', skip = 1)
  %>% rename(Timestamp = Structure)
  %>% mutate(Timestamp = mdy_hm(Timestamp))
  %>% mutate(across(-Timestamp, as.numeric))
)

stopifnot(length(header) == length(names(sm.ca)) - 1)

structure_sm.ca <- tibble(Structure=names(select(sm.ca, where(is.numeric))),
                       Treatment=header)

sm.ca <- (sm.ca
       %>% pivot_longer(where(is.numeric), names_to = "Structure", values_to = "Moisture")
       %>% merge(structure_sm.ca, by="Structure")
       %>% mutate(Date = floor_date(Timestamp, unit = "day"))
       %>% select(Date, Treatment, Moisture)
       %>% group_by(Treatment, Date)
       %>% summarize(Moisture = mean(Moisture, na.rm = TRUE))
       %>% mutate(Treatment=factor(Treatment,
                                   levels=c('50%', '125%'),
                                   ordered=T))
       %>% arrange(Treatment, Date)
)

## produce CA-only figure:
sm.ca.fig <- (
  ggplot(sm.ca, aes(
    x = Date,
    y = Moisture,
    color = Treatment
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
  
)


## WY (all values in m3/m3 VWC)
sm.wy21 <- read.csv("data/soil_moisture/soil_moisture_wy/Restoration_soilmoisture_2021.csv")
sm.wy21 <- sm.wy21 %>% mutate(across(-1, as.numeric))
sm.wy22 <- read.csv("data/soil_moisture/soil_moisture_wy/Restoration_soilmoisture_2022.csv")
sm.wy22 <- sm.wy22 %>% mutate(across(-1, as.numeric))
sm.wy23 <- read.csv("data/soil_moisture/soil_moisture_wy/Restoration_soilmoisture_2023(toAug).csv")
sm.wy23 <- sm.wy23 %>% mutate(across(-1, as.numeric))
#sm.wy24 <- read.csv("data/soil_moisture_wy/Restoration_soilmoisture_2024.csv")
sm.wy <- bind_rows(sm.wy21,sm.wy22,sm.wy23) #combine all annual soil moisture data
sm.wy <- sm.wy %>% mutate(Timestamp = mdy_hm(Measurement.Time))

#pivot longer and combine timestamps into daily soil moisture value
sm.wy.long <- sm.wy %>% pivot_longer(where(is.numeric), names_to = "probenumber", values_to = "Moisture")
sm.wy.long$Block <- gsub("X(\\d+)_Port\\d+", "\\1", sm.wy.long$probenumber)  # Extract block number
sm.wy.long$Port <- gsub(".*_Port([0-9]+)", "\\1", sm.wy.long$probenumber)   # Extract port number
# drought blocks
droughtblocks <- data.frame(Block= c(20,22,27,30,33,42,49,50,51,58),
                            trt = c("drt","drt","cntl","drt","cntl","cntl","cntl","drt","cntl","drt"))
sm.wy.long <- merge(sm.wy.long, droughtblocks)

sm.wy.long <- (sm.wy.long %>% mutate(Date = floor_date(Timestamp, unit = "day"))
         %>% select(Date, trt, Moisture,Block)
         %>% group_by(trt, Date)
         %>% summarize(Moisture = mean(Moisture, na.rm = TRUE))
         %>% arrange(trt, Date)
  )

## make WY-only figure equivalent to figure on 
## CA-only script (make_soil_moistire_figure.R)
sm.wy.fig <- ggplot(sm.wy.long, aes(
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

# # combine if desired. Not cleaned to be pretty and legible yet. 
# ggarrange(gg,smwy.fig, nrow=2)


#### create pretty soil moisture figure with both sites 
## clean up to merge, combine sites, create figure:
sm.wy.long$site <- "Wyoming" #add column to define site
smwy <- sm.wy.long %>% select(trt,Date,Moisture,site) #keep only relevant variables
names(smwy)[names(smwy) == 'trt'] <- 'Treatment' #rename trt variable to match CA
smwy$Treatment <- as.factor(smwy$Treatment) #make factor

sm.ca$Treatment <- factor(sm.ca$Treatment , ordered = FALSE ) #order CA trt factor (to match order of WY when plotting)
sm.ca$site <- "California" #add column to define site

#bring WY and CA data together
allsm <- bind_rows(smwy,sm.ca)

## basic plot (not using)
# ggplot(allsm, (aes(y=Moisture, x=Date, col=Treatment)))+
#   geom_line()+
#   facet_wrap(~site)+
#   scale_color_manual(values = c("blue","skyblue","skyblue","blue"))

## Create pretty plot with correct colors, axes, and shade growing season
## separate out growing seasons in their own dataframes and re-combine
#Create a data.frame for the shaded growing season WY (May - September)
wyshade_data <- data.frame(
  xmin = as.POSIXct(c("2021-05-01", "2022-05-01", "2023-05-01"), tz = "UTC"),
  xmax = as.POSIXct(c("2021-09-01", "2022-09-01", "2023-09-01"), tz = "UTC"),
  ymin = -Inf,  # Extend shading to the bottom of the plot
  ymax = Inf    # Extend shading to the top of the plot
)
#create sequence to use for x axis breaks on both plots
x_breaks <- seq(as.POSIXct("2021-01-01"), as.POSIXct("2023-09-01"), by = "6 months")  # Example: every 6 months

# make WY plot
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
        legend.spacing.y = unit(0.2, "cm"),
        axis.ticks.x = element_line(size=.5))

#Create a data.frame for the shaded growing season CA (January - May)
cashade_data <- data.frame(
  xmin = as.POSIXct(c("2021-01-01", "2022-01-01", "2023-01-01"), tz = "UTC"),
  xmax = as.POSIXct(c("2021-05-01", "2022-05-01", "2023-05-01"), tz = "UTC"),
  ymin = -Inf,  # Extend shading to the bottom of the plot
  ymax = Inf    # Extend shading to the top of the plot
)
sm.ca$Treatment <- ordered(sm.ca$Treatment, levels = c("125%","50%")) #order treatment levels correctly for plot

# make WY plot
casmplot <- ggplot(sm.ca, (aes(y=Moisture, x=Date, col=Treatment)))+
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
        legend.spacing.y = unit(0.2, "cm"),
        axis.ticks.x = element_line(size=.5))

# smwy$Treatment <- relevel(smwy$Treatment, ref = "cntl")
# df$Treatment <- relevel(df$Treatment, ref = "125%")

## make combined figure with axes label
smcombo<- ggarrange(wysmplot,casmplot, nrow=2)
smcombo<-annotate_figure(smcombo, left=text_grob(expression(
  Soil ~ volumetric ~ water ~ content ~ (m ^ 3 ~ m ^ -3)
),rot=90))
smcombo #view

# save as tiff
tiff("figures/sm_figure.tiff", res=400, height = 4,width =7, "in",compression = "lzw")
smcombo
dev.off()

#### Run analysis to see average difference in soil moisture between 
#### precipitation treatments at each site. 
#### ONLY growing season differences are considered.
## WY
## subset data for only growing season
wygrowing <- smwy %>% filter(Date >= as.POSIXct("2021-05-01") & Date <= as.POSIXct("2021-09-01") |
                               Date >= as.POSIXct("2022-05-01") & Date <= as.POSIXct("2022-09-01")|
                               Date >= as.POSIXct("2023-05-01") & Date <= as.POSIXct("2023-09-01"))
## t test for mean differences
t.test(wygrowing$Moisture~wygrowing$Treatment)
## calculate precent reduction to report
((.1838059-.2019800)/.2018900)*100 # percent reduction WY (annual 9%)
((0.2147175-0.2367141)/0.2367141)*100 # percent reduction WY (growing season 9%)

## CA
## subset data for only growing season
cagrowing <- sm.ca %>% filter(Date >= as.POSIXct("2021-01-01") & Date <= as.POSIXct("2021-05-01") |
                               Date >= as.POSIXct("2022-01-01") & Date <= as.POSIXct("2022-05-01")|
                               Date >= as.POSIXct("2023-01-01") & Date <= as.POSIXct("2023-05-01"))
## t test for mean differences
t.test(cagrowing$Moisture~cagrowing$Treatment)
## calculate precent reduction to report
((.1596717-.1819071)/.1819071)*100 # percent reduction CA (annual 12%)
((0.1814996-0.2242060)/0.2242060)*100 # percent reduction WY (growing season 19%)
