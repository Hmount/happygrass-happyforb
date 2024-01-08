#### CWM analysis, 
#### how to CWM's correlate with native cover, invasive cover, FD, and ecosystem services 
## CWM as a predictor of native cover, invasive cover, and ecosystem services

## packages
library(tidyverse)

## load in data, clean and modify columns
# WY
comp.wy <- read.csv("data/comp_wy.csv") # Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
#correct factor columns
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
cwm.wy <- read.csv("data/cwm_wy.csv")# Wyoming CWM data
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities
# combine to master df (remove spp columns for now)
allwy <- merge(comp.wy[,-c(11:66)],cwm.wy, by=c("year","trt","block","subplot"))

# CA
comp.ca <- read.csv("data/comp_ca.csv") #read in California comp data
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
comp.ca$trt <- tolower(comp.ca$trt) #make these lower to match cwm dataframe
comp.ca <- comp.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
cwm.ca <- read.csv("data/cwm_ca.csv")# California CWM data
cwm.ca <- cwm.ca %>% filter(year != "0") #remove predicted/target communities

# combine to master df (remove spp columns for now)
allca <- merge(comp.ca[,-c(18:56)],cwm.ca, by=c("year","trt","block","water"))

