#### CWM summary stats and basic analysis, 
#### how do CWM's maintain/change through time? (dissimilarity analysis)

## packages
library(tidyverse)
#more???

## load in CWM data
cwm.wy <- read.csv("data/cwm_wy.csv")
cwm.wy$year <- as.factor(cwm.wy$year)
#make new sequence column
cwm.wy <- cwm.wy %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.wy$yrorder <- as.numeric(cwm.wy$yrorder)
#add plot ID column (but give NA to seedling/predicted communities)
cwm.wy <- cwm.wy %>% 
  mutate(plot = ifelse(!is.na(subplot), paste(block, trt, subplot, sep = "."), NA))

# CSV of species-trait combinations (for OG 25)
traits.wy <- read.csv("data/mixedgrass.csv", header=TRUE, row.names=1)
traits.wy <- traits.wy[traits.wy$use==1,] # subset use=1
traits.wy$PLSg.m2.mono <- traits.wy$PLSlbs.acre.mono * (453.59237 / 4046.8564224) #convert lb/acre to g/m2

#scale traits
traits.wy$srl = scale(log(traits.wy$srl))
traits.wy$ldmc = scale(log(traits.wy$ldmc))
traits.wy$leafn = scale(log(traits.wy$leafn))
traits.wy$lop = scale(traits.wy$lop)
traits.wy$rootdiam = scale(log(traits.wy$rootdiam))
traits.wy$sla = scale(log(traits.wy$sla))
traits.wy$rdmc = scale(log(traits.wy$rdmc))

# CA
cwm.ca <- read.csv("data/cwm_ca.csv")
cwm.ca$trt <- as.factor(cwm.ca$trt)
cwm.ca$year <- as.factor(cwm.ca$year)
#add plot ID column (but give NA to seedling/predicted communities)
cwm.ca <- cwm.ca %>% 
  mutate(plot = paste(block, trt, sep = "."), NA)

#make new sequence column
cwm.ca <- cwm.ca %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.ca$yrorder <- as.numeric(cwm.ca$yrorder)


library(emmeans)
#### CWM's change overtime
cwm.wy2 <- cwm.wy %>% filter(year!="0") # subset predicted/target data 
## Did CWM's differ significantly across years? (anova's) 
## ~ yes, mostly groups differ signifigantly
m1 <- aov(leafn ~ year*trt, data=cwm.wy2)
TukeyHSD(m1)
m2 <- aov(srl ~ year*trt, data=cwm.wy2)
TukeyHSD(m2)
m3 <- aov(ldmc ~ year*trt, data=cwm.wy2)
TukeyHSD(m3)
m4 <- aov(lop ~ year*trt, data=cwm.wy2)
TukeyHSD(m4)

cwm.ca2 <- cwm.ca %>% filter(year!="0") # subset predicted/target data 
cwm.ca2$trt <- as.factor(cwm.ca2$trt)
cwm.ca2$year <- as.factor(cwm.ca2$year)
## Did CWM's differ significantly across years? (anova's) 
## ~ yes, mostly groups differ signifigantly
m1c <- aov(N ~ year*trt, data=cwm.ca2)
TukeyHSD(m1c)
m2c <- aov(SRL ~ year*trt, data=cwm.ca2)
TukeyHSD(m2c)
m3c <- aov(LMA ~ year*trt, data=cwm.ca2)
TukeyHSD(m3c)
m4c <- aov(seed.mass ~ year*trt, data=cwm.ca2)
TukeyHSD(m4c)

## How did CWM's change overtime? (line plots)
# change in X trait CWM over 3 years (all together) 
# for some traits there may be general trends 
ggplot(cwm.wy2, aes(x=yrorder, y=rootdiam, group=plot, color=trt))+
  geom_point()+
  geom_line()
ggplot(cwm.ca2, aes(x=year, y=SRL, group=plot, color=trt))+
  geom_point()+
  geom_line()

# change in X trait CWM over 3 years separated by seeding trt 
# looks like seeding trt may influence how much CWMs change overtime
ggplot(cwm.wy2, aes(x=yrorder, y=leafn, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_wrap(~trt)

# change in X trait CWM over 3 years separated by seeding trt AND drought trt
# drought also seems to influence how much CWM's changed 
ggplot(cwm.wy2, aes(x=yrorder, y=leafn, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_grid(trt~drought)

## did plots fluxuate a lot? Were certain communities more stable? (how to test?)
## Did certain traits fluxuate more or less? (how to test?)

#### CWM dissimilarity
## NMDS on community trt differences preformed in "CWM_trait_calculations.R"
## how similar are plots to target annually? 
## how different are plots year to year?

# absolute differences between plot CWM and targets
# this is the only way to compare all communities to a single target (other dissimilarity
# could still be calculated for composition data as well as similarity between all plots
# based on trt (also NMDS on other script), or within trt by block or drought)
# how dissimilar are drought and control blocks of the same trt? (effect of drought on CWM) 

# define trait targets
quantile(traits.wy$leafn,.25) #IR
# -0.565257 
quantile(traits.wy$srl,.7557) #IR
# 0.6536232 
quantile(traits.wy$ldmc,.75) # DT
# 0.6766338 
quantile(traits.wy$lop,.25) #DT
# -0.7420366 

## load in data, clean and modify columns (check this)
# WY
comp.wy <- read.csv("data/hpg_total.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)

comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data

comp.wy.wide$plot.tveg <- comp.wy.wide$plot.tveg/100  # make native live veg % a proportion to match CA data

# merge with CWM data
cwm.wy <- read.csv("data/cwm_wy") #cwm data
alldat.wy <- merge(comp.wy.wide, cwm.wy, by = "block") #combine



# Determine how far off from the targets we were for EVERY treatment (refined in next steps)
cwm.wy$dist_leafn <- abs(cwm.wy$leafn - quantile(traits.wy$leafn,.25))
cwm.wy$dist_srl <- abs(cwm.wy$srl - quantile(traits.wy$srl,.7557))
cwm.wy$dist_ldmc <- abs(cwm.wy$ldmc - quantile(traits.wy$ldmc,.75))
cwm.wy$dist_lop <- abs(cwm.wy$lop - quantile(traits.wy$lop,.25))

# Calculate how far we were from the dt, ir, and fd objective. fd objective is simply highest diversity after normalization.
cwm.wy$dist_dt <- rowSums(cbind(cwm.wy$dist_ldmc,cwm.wy$dist_lop))
cwm.wy$dist_ir <- rowSums(cbind(cwm.wy$dist_leafn,cwm.wy$dist_srl))
cwm.wy_roaq <- read.csv("data/cwm_raoq_wy.csv")
cwm.wy <- merge(cwm.wy, cwm.wy_roaq[,c(7:12)], all=T)
# cwm.wy_roaq <- cwm.wy_roaq
# cwm.wy$raoq <- cwm.wy_roaq$full

ggplot(cwm.wy, aes(y=sqrt(nativecov), x=srl, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~year)

ggplot(cwm.wy, aes(y=sqrt(nativecov), x=dist_dt, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~year)

allwy$dist_leafn <- abs(allwy$leafn - quantile(traits.wy$leafn,.25))
allwy$dist_srl <- abs(allwy$srl - quantile(traits.wy$srl,.7557))
allwy$dist_ldmc <- abs(allwy$ldmc - quantile(traits.wy$ldmc,.75))
allwy$dist_lop <- abs(allwy$lop - quantile(traits.wy$lop,.25))

# Calculate how far we were from the dt, ir, and fd objective. fd objective is simply highest diversity after normalization.
allwy$dist_dt <- rowSums(cbind(allwy$dist_ldmc,allwy$dist_lop))
allwy$dist_ir <- rowSums(cbind(allwy$dist_leafn,allwy$dist_srl))
allwy_roaq <- read.csv("data/cwm_raoq_wy.csv")
allwy <- merge(allwy, allwy_roaq[,c(7:12)], all=T)

justdt <- allwy %>% filter(trt=="dt")
ggplot(justdt, aes(y=sqrt(nativecov), x=dist_dt, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought.y~year)

justir <- allwy %>% filter(trt=="ir")
ggplot(justir, aes(y=sqrt(nativecov), x=dist_ir, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought.y~year)

ggplot(allwy, aes(y=sqrt(nativecov), x=full, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought.y~year)

# # euclidean distances (dissimilarity) between plot CWM and targets
# # could do this with trt group means, but why? previous method is only applicable to 
# # trait values (although compositional dissimilarity could be examined)
# #subset communities by trt
# allsub <- allwy %>% filter(trt=="ir")
# #calculate distance from target for relvent traits
# #ir.cwm.dist <- dist(allsub$)
# #recombine
# 
# #dissimilarity between years 
# library(vegan)
# # find dissimilrity to tareget for each plot each year
# # in loop for each trait for each plot 
# # subset data for one year at a time 
# # make into matrixes
# # calculate dist from target
# # save in dataframe
# # Calculate dissimilarity for each year






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

