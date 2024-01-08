#### CWM summary stats and basic analysis, 
#### how do CWM's maintain/change through time? (dissimilarity analysis)

## packages
library(tidyverse)
#more???

## load in CWM data
cwm.wy <- read.csv("data/cwm_wy.csv")
#make new sequence column
cwm.wy <- cwm.wy %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.wy$yrorder <- as.numeric(cwm.wy$yrorder)
#add plot ID column (but give NA to seedling/predicted communities)
cwm.wy <- cwm.wy %>% 
  mutate(plot = ifelse(!is.na(subplot), paste(block, trt, subplot, sep = "."), NA))

cwm.ca <- read.csv("data/cwm_ca.csv")
#make new sequence column
cwm.ca <- cwm.ca %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.ca$yrorder <- as.numeric(cwm.ca$yrorder)


#### How did CWM's change overtime? 
## did trt's differ annually? 
## Did certain traits fluxuate more or less? (how to test?)
## did plots fluxuate a lot? (how to test?)
ggplot(test, aes(x=yrorder, y=leafn), color=plot)+
  geom_point()+
  geom_line()
  #geom_smooth(method = "lm")
test <- cwm.wy %>% filter(year!="0") 
test[,c(11:18)] <- as.factor(test[,c(11:18)])
ggplot(test, aes(x=yrorder, y=sla, group=plot), color=trt)+
  geom_point()+
  geom_line()#+
facet_wrap(~yrorder)
test <- cwm.wy %>% filter(subplot =="s") 
test$plot <- as.factor(test$plot)

plot(cwm.wy$sla~cwm.wy$yrorder)
test[,c(11:18)] <- lapply(test[,c(11:18)], as.factor)

ggplot(test, aes(x=yrorder, y=sla, group=plot), color=trt)+
  geom_point()+
  geom_line()+
  scale_color_identity(name="trt",
                       labels=c("dt","ir","rand","fd"), 
                       guide="legend")

ggplot(test, aes(x=yrorder, y=leafn, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_wrap(~trt)

ggplot(test, aes(x=yrorder, y=rootdiam, group=plot, color=trt))+
  geom_point()+
  geom_line()

ggplot(test, aes(x=yrorder, y=leafn, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_grid(trt~drought)


#### CWM dissimilarity
## how similar are plots to target annually? 
## how different are plots year to year?

#dissimilarity between years 
library(vegan)
# Calculate dissimilarity for each year
dissimilarity_list <- lapply(1:3, function(i) {
  target <- cwm_p.wy[, 2:9]  # Assuming columns 2:8 are the traits in preds
  actual <- get(paste0("cwm2", i))[, 2:9]  # Assuming columns 2:8 are the traits in cwm21, cwm22, cwm23
  calculate_bray_curtis(target, actual)
})
# Assuming preds, cwm21, cwm22, cwm23 are your dataframes

# Define a function to calculate Bray-Curtis dissimilarity
calculate_bray_curtis <- function(target, actual) {
  vegdist(data.frame(target, actual), method = "bray")
}

# Calculate dissimilarity for each year
dissimilarity_list <- lapply(1:3, function(i) {
  target <- preds[, 2:8]  # Assuming columns 2:8 are the traits in preds
  actual <- get(paste0("cwm2", i))[, 2:8]  # Assuming columns 2:8 are the traits in cwm21, cwm22, cwm23
  calculate_bray_curtis(target, actual)
})
?vegdist
# Bind the dissimilarity matrices into a single dataframe
dissimilarity_df <- do.call(cbind, dissimilarity_list)

# Attach PlotID (assuming it's a column in preds) to the dissimilarity matrix
dissimilarity_df <- cbind(PlotID = preds$PlotID, dissimilarity_df)


ggplot(cwm.wy, aes(x=trt, y=leafn), color=yrorder)+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(cwm.wy, aes(x=block, y=leafn), color=yrorder)+
  geom_point()+
  geom_smooth(method = "lm")
ggplot(cwm.wy, aes(x=yrorder, y=leafn), color=trt)+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~trt)
ggplot(cwm.wy, aes(x=yrorder, y=leafn), color=block)+
  geom_point()+
  geom_smooth(method = "lm")

## load in data, clean and modify columns (check this)
# WY
comp.wy <- read.csv("data/hpg_total.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
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
rm(wy.invaded.23)# remove invasion data
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
comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data

# calculate cover native per PLOT (averaged per subplot)
forplotcover <- comp.wy %>% group_by(year,block,trt,species) %>% 
  mutate(species.plotmean = mean(cover)) #mean per species per PLOT
forplotcover2 <- forplotcover %>% filter(species!="BG"&
                                           species!="Litter"&
                                           native == "N" &
                                           subplot == "n") %>% #only native live veg
  group_by(year,block,trt) %>% 
  summarize(plot.tveg = sum(species.plotmean, na.rm=T)) #summarize total live veg per PLOT
comp.wy.wide <- merge(comp.wy.wide,forplotcover2, by=c("year","block","trt"),all.x = T)
comp.wy.wide$plot.tveg <- comp.wy.wide$plot.tveg/100  # make native live veg % a proportion to match CA data

# merge with CWM data
cwm.wy <- read.csv("data/cwm_wy") #cwm data
alldat.wy <- merge(comp.wy.wide, cwm.wy, by = "block") #combine


# CA
comp.ca <- read.csv("data/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)] #remove empty columns
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)

# merge with CWM data
cwm.ca <- read.csv("data/cwm_ca") #cwm data
alldat.ca <- merge(comp.ca.wide, cwm.ca, by = "trt.b") #combine

