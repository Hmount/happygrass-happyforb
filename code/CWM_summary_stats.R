#### CWM summary stats and basic analysis, 
#### How do CWM's maintain/change through time? 
#### Are CWM's effected by drought and invasion treatments (and seeding trt?)
#### Dissimilarity analysis

## packages
library(tidyverse)

#### load in CWM data and raw trait matricies ####
## WY
cwm.wy <- read.csv("data/cwm_wy.csv") #CWM data
cwm.wy$year <- as.factor(cwm.wy$year)
cwm.wy <- cwm.wy %>% mutate(yrorder = ifelse(year=="2021","1", #make new sequence column
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.wy$yrorder <- as.numeric(cwm.wy$yrorder)
cwm.wy <- cwm.wy %>% #add plot ID column (but give NA to target/predicted communities)
  mutate(plot = ifelse(!is.na(subplot), paste(block, trt, subplot, sep = "."), NA))
#trait data
traits.wy <- read.csv("data/mixedgrass.csv", header=TRUE, row.names=1) # CSV of species-trait combinations (for OG 25)
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

## CA
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
#traits
traits.ca <- read.csv("data/annualgrass.csv", header=TRUE, row.names=1) # CSV of species-trait combinations (for OG 25)
traits.ca <- traits.ca %>% select(-c(Asat,WUE,RLD,RTD)) #remove unused traits
traits.ca$graminoid <- c(0,0,0,0,1,0,1,0,0,0, #10   # add graminoid
                         1,1,0,1,0,1,1,1,0,0, #20
                         0,0,0,0,1,1,0,0,0,1, #30
                         1,0,0,1,1,0,1,1,0,1,1) #41
#scale traits
traits.ca$SRL = scale(log(traits.ca$SRL))
traits.ca$LMA = scale(log(traits.ca$LMA))
traits.ca$N = scale(log(traits.ca$N))
traits.ca$seed.mass = scale(traits.ca$seed.mass)
traits.ca$Rdiam = scale(log(traits.ca$Rdiam))
#traits.ca$sla = scale(log(traits.ca$sla))
traits.ca$RMF = scale(log(traits.ca$RMF))
#filter for only species used (OG 25)
traits.ca <- traits.ca %>% 
  filter(Code != "AVBA") %>% 
  filter(Code != "BRNI") %>%
  filter(Code != "LUAL") %>%
  filter(Code != "MASA") %>% 
  filter(Code != "PEHE") %>% 
  filter(Code != "SACO") %>%
  filter(Code != "FEPE") %>%  # this one appears in communities, just not OG preds
  filter(Code != "FEMY") %>%  # this one appears in communities, just not OG preds
  filter(Code != "BRMA")      # this one appears in communities, just not OG preds

leafn.21.mod <- aov(leafn~trt, cwm21.wy)



#### CWM's change overtime ####
## WY
cwm.wy2 <- cwm.wy %>% filter(year!="0") # subset predicted/target data 
## Did CWM's differ significantly across years? (anova's) 
## ~ yes, mostly groups differ signifigantly
leafn.m1 <- aov(leafn ~ year*trt, data=cwm.wy2)
summary(leafn.m1) #trt + year + int sig
TukeyHSD(leafn.m1)
srl.m1 <- aov(srl ~ year*trt, data=cwm.wy2)
summary(srl.m1) #trt + year + int sig
TukeyHSD(srl.m1)
veg.m1 <- aov(veg ~ year*trt, data=cwm.wy2)
summary(veg.m1) #trt + year + int sig
TukeyHSD(veg.m1)
ldmc.m1 <- aov(ldmc ~ year*trt, data=cwm.wy2)
summary(ldmc.m1) #trt + year + int sig
TukeyHSD(ldmc.m1)
lop.m1 <- aov(lop ~ year*trt, data=cwm.wy2)
summary(lop.m1) #trt + year + int sig
TukeyHSD(lop.m1)
rd.m1 <- aov(rootdiam ~ year*trt, data=cwm.wy2)
summary(rd.m1) #trt + year + int sig
TukeyHSD(rd.m1)

fd.m1 <- aov(rootdiam ~ year*trt, data=cwm.wy2)
summary(rd.m1) #trt + year + int sig
TukeyHSD(rd.m1)

## CA
cwm.ca2 <- cwm.ca %>% filter(year!="0") # subset predicted/target data 
cwm.ca2$trt <- as.factor(cwm.ca2$trt)
cwm.ca2$year <- as.factor(cwm.ca2$year)
## Did CWM's differ significantly across years? (anova's) 
## ~ yes, mostly groups differ signifigantly
N.m1 <- aov(N ~ year*trt, data=cwm.ca2)
summary(N.m1) #trt + year sig
TukeyHSD(N.m1)
SRL.m1 <- aov(SRL ~ year*trt, data=cwm.ca2)
summary(SRL.m1) #trt + year sig
TukeyHSD(SRL.m1)
RMF.m1 <- aov(RMF ~ year*trt, data=cwm.ca2)
summary(RMF.m1) #not sig
TukeyHSD(RMF.m1)
LMA.m1 <- aov(LMA ~ year*trt, data=cwm.ca2)
summary(LMA.m1) #trt + year + interaction sig
TukeyHSD(LMA.m1)
SM.m1 <- aov(seed.mass ~ year*trt, data=cwm.ca2)
summary(SM.m1) #not sig
TukeyHSD(SM.m1)
RD.m1 <- aov(Rdiam ~ year*trt, data=cwm.ca2)
summary(RD.m1) #not sig
TukeyHSD(RD.m1)

## How did CWM's change overtime? (line plots)
## change in X trait CWM over 3 years (all together) 
# for some traits there may be general trends 
ggplot(cwm.wy2, aes(x=yrorder, y=rootdiam, group=plot, color=trt))+
  geom_point()+
  geom_line()
ggplot(cwm.ca2, aes(x=year, y=SRL, group=plot, color=trt))+
  geom_point()+
  geom_line()

## change in X trait CWM over 3 years separated by seeding trt 
ggplot(cwm.wy2, aes(x=yrorder, y=leafn, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_wrap(~trt)
#looks like seeding trt may influence how much CWMs change overtime
ggplot(cwm.ca2, aes(x=yrorder, y=N, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_wrap(~trt)
#looks like similar traends overtime regardless of seeding trt

## change in X trait CWM over 3 years separated by seeding trt AND drought trt
## these (i.e. effect of drought) are not statistically modeled (yet)
ggplot(cwm.wy2, aes(x=yrorder, y=leafn, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_grid(trt~drought)
# drought also seems to influence how much CWM's changed 
ggplot(cwm.ca2, aes(x=yrorder, y=N, group=plot, color=trt))+
  geom_point()+
  geom_line()+
  facet_grid(trt~water)
# drought may have had effects, but not by trt 

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
comp.wy <- read.csv("data/comp_wy.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy$nativecov <- comp.wy$nativecov/100  # make native live veg % a proportion to match CA data
comp.wy$plot.tveg <- comp.wy$plot.tveg/100  # make native live veg % a proportion to match CA data
# merge with CWM data
cwm.wy$subplot <- as.factor(cwm.wy$subplot)
cwm.wy$drought <- as.factor(cwm.wy$drought)
cwm.wy$year<-as.factor(cwm.wy$year)
cwm.wy$trt <- as.factor(cwm.wy$trt)
cwm.wy$block <- as.factor(cwm.wy$block)
#alldat.wy <- merge(comp.wy, cwm.wy, by = c("year", "block", "trt", "subplot")) #combine

alldat.wy <- merge(comp.wy, cwm.wy, by = c("year", "block", "trt", "subplot","plot")) #combine

#distdata
#merge

ggplot(cwm.wy, aes(y=sqrt(nativecov), x=srl, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~year)

ggplot(cwm.wy, aes(y=sqrt(nativecov), x=dist_dt, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought~year)


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

cwm.wy <- read.csv("data/cwm_wy.csv")# 
cwm.wy2 <- read.csv("data/cwm_wy(plot).csv")#




# CA
comp.ca <- read.csv("data/comp_ca.csv") #read in California comp data
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
# merge with CWM data
cwm.ca$water <- as.factor(cwm.ca$water)
alldat.ca <- merge(comp.ca, cwm.ca, by = "trt.b") #combine

