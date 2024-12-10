library(tidyverse)

header <- (scan("data/soil_moisture/soil_moisture_ca/three_year_soil_moisture.csv",
                nlines = 1,
                sep=",",
                what = character())
           %>% tail(-1))

df <- (
  readr::read_csv('data/soil_moisture/soil_moisture_ca/three_year_soil_moisture.csv', skip = 1)
  %>% rename(Timestamp = Structure)
  %>% mutate(Timestamp = mdy_hm(Timestamp))
  %>% mutate(across(-Timestamp, as.numeric))
)

stopifnot(length(header) == length(names(df)) - 1)

structure_df <- tibble(Structure=names(select(df, where(is.numeric))),
                      Treatment=header)

df <- (df
  %>% pivot_longer(where(is.numeric), names_to = "Structure", values_to = "Moisture")
  %>% merge(structure_df, by="Structure")
  %>% mutate(Date = floor_date(Timestamp, unit = "day"))
  %>% select(Date, Treatment, Moisture)
  %>% group_by(Treatment, Date)
  %>% summarize(Moisture = mean(Moisture, na.rm = TRUE))
  %>% mutate(Treatment=factor(Treatment,
                              levels=c('50%', '125%'),
                              ordered=T))
  %>% arrange(Treatment, Date)
)

gg <- (
  ggplot(df, aes(
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

start <- as.POSIXct(strptime("2021-01-01 00:00", format = "%Y-%m-%d %H:%M"))
specific.breaks <- c(start)
for (ix in 1:5) {
  prev <- tail(specific.breaks, 1)
  specific.breaks <- c(specific.breaks, prev + months(6))
}
show(gg + scale_x_datetime(labels = scales::date_format("%Y-%m"), 
                           breaks = specific.breaks,
                           #oob=scales::oob_keep()
                           limits = lims 
                                  #expand = c(0, 0)
                           ))


ggsave(
  'soil-moisture-levels-over-time.pdf',
  plot = gg,
  height = 8,
  width = 12,
  units = 'in'
)

ggsave(
  'soil-moisture-levels-over-time.tiff',
  plot = gg,
  height = 8,
  width = 12,
  units = 'in',
  dpi = 600
)
