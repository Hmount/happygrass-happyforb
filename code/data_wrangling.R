#### Data wrangling
#### take raw data (or compiled data from Jen/CA) and wrangling 

####
#### CA
## load in data, clean and modify columns
comp.ca <- read.csv("data/raw_cover/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)] #remove empty columns
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)

# save data for all other analyses (changing columns to factors is not preserved so this needs to be repeated)
write.csv(comp.ca, "data/comp_ca.csv", row.names = F)


####
#### WY (subplot level)
## load in data, clean and modify columns
comp.wy <- read.csv("data/raw_cover/hpg_total.csv") #Wyoming species comp data 
#comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
#correct factor columns
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy <- comp.wy %>% select(-c(sub.tcov,sub.tveg)) #remove totals calculated in excel so all math is recorded here

# add invaded locations as 0/1 variable (2023 only)
wy.invaded.23 <- read.csv("data/invasion_loc23.csv") # load data
wy.invaded.23$invaded <-tolower(wy.invaded.23$invaded)
comp.wy <- comp.wy %>%
  mutate(invaded = ifelse(year == 2023, as.integer(paste(block, trt, subplot) %in% paste(wy.invaded.23$block, wy.invaded.23$trt, wy.invaded.23$invaded)), NA_integer_))
#rm(wy.invaded.23)# remove invasion data

# make unique plot variable
comp.wy <- comp.wy %>% unite(plot, c(block, trt, subplot), sep = ".", remove=F) 

## Summaries at the subplot level:
# calculate cover per subplot and add to data
fortotalcover <- comp.wy %>% filter(species!="BG"&
                                       species!="Litter") %>% #only native live veg (this is unecessary as I already seperated out BG)
  group_by(year,block,trt,subplot) %>% 
  summarize(totalcov = sum(cover, na.rm=T)) #summarize total live veg per subplot
comp.wy <- merge(comp.wy,fortotalcover, all.x = T)

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
                                                         "nativecov","totalcov","BG", "Litter",
                                                         "plot", "invaded"), 
                                             names_from = "species", 
                                             values_from = "cover")
comp.wy.wide$totalcov <- comp.wy.wide$totalcov/100  # make native live veg % a proportion to match CA data
comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data

# save data for all other analyses
write.csv(comp.wy.wide, "data/comp_wy.csv", row.names = F)


####
#### WY (again for PLOT level) 
#### invaded is not relevent at plot level so not added to this df
comp.wy.plot <- read.csv("data/raw_cover/hpg_total.csv") #Wyoming species comp data 
#comp.wy.plot <- comp.wy.plot %>% filter(year != "2020") #remove pre-treatment data
#correct factor columns
comp.wy.plot$subplot <- as.factor(comp.wy.plot$subplot)
comp.wy.plot$drought <- as.factor(comp.wy.plot$drought)
comp.wy.plot$year<-as.factor(comp.wy.plot$year)
comp.wy.plot$trt <- as.factor(comp.wy.plot$trt)
comp.wy.plot$block <- as.factor(comp.wy.plot$block)
comp.wy.plot <- comp.wy.plot %>% select(-c(sub.tcov,sub.tveg)) #remove totals calculated in excel so all math is recorded here

# make unique plot variable
comp.wy.plot <- comp.wy.plot %>% unite(plot, c(block, trt), sep = ".", remove=F) 

# calculate total cover per PLOT (averaged per subplot)
fortotplotcover <- comp.wy.plot %>% group_by(year,block,trt,species) %>% 
  mutate(tot.plotmean = mean(cover)) #mean per species per PLOT (using mutate to perserve columns, but averages appear twice)

#first summarize litter and BG
forbg <- comp.wy.plot %>% group_by(year,block,trt,drought) %>% 
  summarize(bg.plotmean = mean(BG))
forlit <- comp.wy.plot %>% group_by(year,block,trt,drought) %>% 
  summarize(lit.plotmean = mean(Litter))

#now summarize plants
fortotplotcover <- comp.wy.plot %>% group_by(year,block,trt,drought,species) %>% 
  summarize(tot.plotmean = mean(cover))
plotcovers <- fortotplotcover %>% pivot_wider(id_cols = c("year","block","trt","drought"), 
                                                       names_from = "species", 
                                                       values_from = "tot.plotmean")
fortotplotcover2 <- fortotplotcover %>% group_by(year,block,trt) %>% 
  summarize(totalcov.plot = sum(tot.plotmean, na.rm=T)) #summarize total live veg per PLOT
#combine
comp.wy.plot.wide <- merge(fortotplotcover2,plotcovers) 
forplotcover2 <- forplotcover %>% group_by(year,block,trt,drought,species) %>% 
  filter(native == "N") %>% #only native live veg
  group_by(year,block,trt) %>% 
  summarize(nativecov.plot = sum(species.plotmean, na.rm=T)) #summarize total live veg per PLOT
comp.wy.plot.wide <- merge(comp.wy.plot.wide,forplotcover2, all = T) #combine
#more merge
comp.wy.plot.wide <- merge(comp.wy.plot.wide,forbg)
comp.wy.plot.wide <- merge(comp.wy.plot.wide,forlit)

# # make wide for analysis and matching CA data
# comp.wy.plot.wide <- comp.wy.plot.wide %>% select(-c("native"))#,"graminoid")) #columns to drop 
# comp.wy.plot.wide <- comp.wy.plot.wide %>% pivot_wider(id_cols = c("year","block","trt","drought",
#                                                                    "nativecov.plot","totalcov.plot"), 
#                                                        names_from = "species", 
#                                                        values_from = "cover")
comp.wy.plot.wide$totalcov.plot <- comp.wy.plot.wide$totalcov.plot/100  # make native live veg % a proportion to match CA data
comp.wy.plot.wide$nativecov.plot <- comp.wy.plot.wide$nativecov.plot/100  # make native live veg % a proportion to match CA data
comp.wy.plot.wide$bg.plotmean <- comp.wy.plot.wide$bg.plotmean/100  # make native live veg % a proportion to match CA data
comp.wy.plot.wide$lit.plotmean <- comp.wy.plot.wide$lit.plotmean/100  # make native live veg % a proportion to match CA data


#checking in averaging went right (it did)
tot <- forplotcover %>% filter(plot == "1.dt.s") %>% #summed cover of south or north
  group_by(year,block,trt,subplot) %>% 
  summarise(total = sum(cover, na.rm=T))
(62.80+59.7)/2 #add together and divide (checked for each year on one example, all correct)
comp.wy.plot.wide %>% filter(plot=="1.dt.n")

# save data for all other analyses
write.csv(comp.wy.plot.wide, "data/comp_wy_plot.csv", row.names = F)


# # BELOW IS ALL GARBAGE #
# 
# fortotplotcover2 <- fortotplotcover %>% filter(species!="BG"&
#                                            species!="Litter") %>%
#   group_by(year,block,trt) %>% 
#   summarize(totalcov.plot = sum(tot.plotmean, na.rm=T)) #summarize total live veg per PLOT
# comp.wy.plot <- merge(comp.wy.plot,fortotplotcover2, all.x = T)
# 
# # calculate cover native per PLOT (averaged per subplot)
# forplotcover <- comp.wy.plot %>% group_by(year,block,trt,species) %>% 
#   mutate(species.plotmean = mean(cover)) #mean per species per PLOT (using mutate to perserve columns, but averages appear twice)
# forplotcover2 <- forplotcover %>% filter(species!="BG"&
#                                            species!="Litter"&
#                                            native == "N" & #only native live veg
#                                            subplot == "n") %>% #need to use only one subplot since averages appear twice
#   group_by(year,block,trt) %>% 
#   summarize(nativecov.plot = sum(species.plotmean, na.rm=T)) #summarize total live veg per PLOT
# comp.wy.plot <- merge(comp.wy.plot,forplotcover2, all.x = T)
# 
# # make wide for analysis and matching CA data
# comp.wy.plot.wide <- comp.wy.plot %>% select(-c("prob","native","graminoid")) #columns to drop 
# comp.wy.plot.wide <- comp.wy.plot.wide %>% pivot_wider(id_cols = c("year","block","trt","subplot","drought",
#                                                          "nativecov.plot","totalcov.plot","BG", "Litter",
#                                                          "plot"), 
#                                              names_from = "species", 
#                                              values_from = "cover")
# comp.wy.plot.wide$totalcov.plot <- comp.wy.plot.wide$totalcov.plot/100  # make native live veg % a proportion to match CA data
# comp.wy.plot.wide$nativecov.plot <- comp.wy.plot.wide$nativecov.plot/100  # make native live veg % a proportion to match CA data


