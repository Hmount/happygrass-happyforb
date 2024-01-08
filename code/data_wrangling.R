#### Data wrangling
#### take raw data (or compiled data from Jen) and wrangling 

## load in data, clean and modify columns
# WY
comp.wy <- read.csv("data/raw_cover/hpg_total.csv") #Wyoming species comp data 
#comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
#correct factor columns
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)

# add invaded locations as 0/1 variable (2023 only)
wy.invaded.23 <- read.csv("data/invasion_loc23.csv") # load data
wy.invaded.23$invaded <-tolower(wy.invaded.23$invaded)
comp.wy <- comp.wy %>%
  mutate(invaded = ifelse(year == 2023, as.integer(paste(block, trt, subplot) %in% paste(wy.invaded.23$block, wy.invaded.23$trt, wy.invaded.23$invaded)), NA_integer_))
#rm(wy.invaded.23)# remove invasion data

# make unique plot variable
comp.wy <- comp.wy %>% unite(plot, c(block, trt, subplot), sep = ".", remove=F) 

# calculate cover native per subplot and add to data
fornativecover <- comp.wy %>% filter(species!="BG"&
                                       species!="Litter"&
                                       native == "N") %>% #only native live veg
  group_by(year,block,trt,subplot) %>% 
  summarize(nativecov = sum(cover, na.rm=T)) #summarize total live veg per subplot
comp.wy <- merge(comp.wy,fornativecover, all.x = T)

# make wide for analysis and matching CA data
comp.wy.wide <- comp.wy %>% select(-c("prob","native","graminoid")) #columns to drop 
comp.wy.wide <- comp.wy.wide %>% pivot_wider(id_cols = c("year","block","trt","subplot","drought",
                                                         "nativecov","BG", "Litter","plot", "invaded"), 
                                             names_from = "species", 
                                             values_from = "cover")
#comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data

# calculate cover native per PLOT (averaged per subplot)
forplotcover <- comp.wy %>% group_by(year,block,trt,species) %>% 
  mutate(species.plotmean = mean(cover)) #mean per species per PLOT (using mutate to perserve columns, but averages appear twice)
forplotcover2 <- forplotcover %>% filter(species!="BG"&
                                           species!="Litter"&
                                           native == "N" & #only native live veg
                                           subplot == "n") %>% #need to use only one subplot since averages appear twice
  group_by(year,block,trt) %>% 
  summarize(plot.tveg = sum(species.plotmean, na.rm=T)) #summarize total live veg per PLOT
comp.wy.wide <- merge(comp.wy.wide,forplotcover2, by=c("year","block","trt"),all.x = T)
#comp.wy.wide$plot.tveg <- comp.wy.wide$plot.tveg/100  # make native live veg % a proportion to match CA data

# save data for all other analyses
write.csv(comp.wy.wide, "data/comp_wy.csv", row.names = F)

comp.wy <- read.csv("data/comp_wy.csv")

tot <- forplotcover %>% filter(plot == "1.dt.s") %>% 
  group_by(year,block,trt,subplot) %>% 
  summarise(total = sum(cover, na.rm=T))
(62.80+59.7)/2

# CA
comp.ca <- read.csv("data/raw_cover/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)] #remove empty columns
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
#make ca long
#comp.ca.long <- comp.ca %>% pivot_longer(cols = c(18:57), 
#                                         names_to = "species", 
#                                         values_to = "cover")

# save data for all other analyses (changing columns to factors is not preserved so this needs to be repeated)
write.csv(comp.ca, "data/comp_ca.csv", row.names = F)

