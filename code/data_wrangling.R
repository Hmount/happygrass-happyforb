#### Data wrangling
#### take raw data (or compiled data from Jen/CA) and do some tidying to get in
#### order to work with for analyses. All covers are converted to relative 
#### proportions of each species. 

####
#### CA
## load in data, clean and modify columns
comp.ca <- read.csv("data/raw_cover/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)] #remove empty columns
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year <-as.factor(comp.ca$Year)
comp.ca$block <- as.factor(comp.ca$block)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$trt <- tolower(comp.ca$trt) #make lowercase to match

# briefly make long to calculate relative cover by plot of all species
comp.ca.long <- comp.ca %>% pivot_longer(c(18:56), names_to = "species",values_to = "abscov")
# group plots together and find plot totals
comp.ca.long <- comp.ca.long %>% 
  group_by(block, trt, year) %>% 
  mutate(totcov.plot = sum(abscov, na.rm=T))
# find relative proportion per species
comp.ca.long <- comp.ca.long %>% 
  group_by(block, trt, year) %>% 
  mutate(relcov = abscov/totcov.plot)
# #checks
# comp.ca.long %>% filter(plot=="1-DT") %>% summarize(totcov.plot = sum(abscov, na.rm=T))
# comp.ca.long %>% filter(plot=="1-DT") %>% summarize(tot = sum(relcov, na.rm=T))
# return to wide
comp.ca.long <- comp.ca.long %>% 
  pivot_wider(id_cols = c(year,block,trt,water,structure), names_from = "species", values_from = "relcov")    #longer(c(18:56), names_to = "species",values_to = "abscov")
comp.ca.rel <- merge(comp.ca[,c(2:17,58)],comp.ca.long,all.y=T, by=c("year","block","trt","water","structure"))
# a few NaN produced in empty communities (36.fd.22, 36.fd.23, and 47.ir.22)
comp.ca.rel[is.na(comp.ca.rel)] <- NA # made into 0 here, but do not analyze these.
# confirm all communities sum to approximately 1 (except 3 empty plots listed above)
rowSums(comp.ca.rel[,c(18:56)], na.rm=T)

## save data for all other analyses 
# (changing columns to factors is not preserved so this needs to be repeated)
# notice that the new dataframe of relative covers is saved and used in all 
# downstream analyses
write.csv(comp.ca.rel, "data/comp_ca.csv", row.names = F)


####
#### WY (subplot level) 
#### NOT UPDATED W/ RELATIVE COV (!)
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
#### invaded is not relevant at plot level so not added to this df
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

# # calculate total cover per PLOT (averaged per subplot)
# fortotplotcover <- comp.wy.plot %>% group_by(year,block,trt,species) %>% 
#   mutate(tot.plotmean = mean(cover)) #mean per species per PLOT (using mutate to perserve columns, but averages appear twice)

# first, summarize means of litter and BG by plot 
forbg <- comp.wy.plot %>% group_by(year,block,trt,drought) %>% 
  summarize(bg.plotmean = mean(BG))
forbg$bg.plotmean <- forbg$bg.plotmean/100  # make into proportion
forlit <- comp.wy.plot %>% group_by(year,block,trt,drought) %>% 
  summarize(lit.plotmean = mean(Litter))
forlit$lit.plotmean <- forlit$lit.plotmean/100  # make into proportion

# now summarize plants
comp.wy.plot$cover <- comp.wy.plot$cover/100  # make into proportion
fortotplotcover <- comp.wy.plot %>%  
  group_by(year,block,trt,drought,species,native,graminoid) %>% 
  summarize(tot.plotmean = mean(cover)) #summarize mean absolute cover/spp/plot 

# #checking in averaging went right (it did)
# tot <- forplotcover %>% filter(plot == "1.dt.s") %>% #summed cover of south or north
#   group_by(year,block,trt,subplot) %>% 
#   summarise(total = sum(cover, na.rm=T))
# (62.80+59.7)/2 #add together and divide (checked for each year on one example, all correct)
# comp.wy.plot.wide %>% filter(plot=="1.dt.n")

# group plots together and find plot totals
fortotplotcover <- fortotplotcover %>% 
  group_by(block, trt, year) %>% 
  mutate(totcov.plot = sum(tot.plotmean, na.rm=T)) #total absolute cover of all plants per plot
fornativeplotcover <- fortotplotcover %>% group_by(year,block,trt,drought) %>% 
  filter(native == "N") %>% #only native live veg
  group_by(year,block,trt) %>% 
  summarise(nativecov.plot = sum(tot.plotmean, na.rm=T)) #sum relative cover live native veg per PLOT

# find relative proportion per species
fortotplotcover <- fortotplotcover %>% 
  group_by(block, trt, year) %>% 
  mutate(relcov = tot.plotmean/totcov.plot)
# #checks
# fortotplotcover %>% filter(block=="1") %>% filter(trt=="dt") %>% summarize(totcov.plot = sum(tot.plotmean, na.rm=T))
# fortotplotcover %>% filter(block=="1"&trt=="dt") %>% summarize(tot = sum(relcov, na.rm=T))

# make into wide format
plotcovers <- fortotplotcover %>% pivot_wider(id_cols = c("year","block","trt","drought"), 
                                    names_from = "species", 
                                    values_from = "relcov")
# combine
comp.wy.plot.wide <- merge(fortotplotcover[,c(1:4,9)],plotcovers)
comp.wy.plot.wide <- merge(comp.wy.plot.wide, fornativeplotcover, all.x=T)
comp.wy.plot.wide <- merge(comp.wy.plot.wide, forbg, all.x=T)
comp.wy.plot.wide <- merge(comp.wy.plot.wide, forlit, all.x=T)
comp.wy.plot.wide <- distinct(comp.wy.plot.wide) #remove duplicate rows from merge
# ensure all NA or NAN are NA that will sun to 0
comp.wy.plot.wide[is.na(comp.wy.plot.wide)] <- NA # made into 0 here, but do not analyze these.
# confirm all communities sum to approximately 1 (except 3 empty plots listed above)
rowSums(comp.wy.plot.wide[,c(6:61)], na.rm=T)

## save data for all other analyses 
# (changing columns to factors is not preserved so this needs to be repeated)
# notice that the new dataframe of relative covers is saved and used in all 
# downstream analyses
#
# save data for all other analyses
write.csv(comp.wy.plot.wide, "data/comp_wy_plot.csv", row.names = F)
