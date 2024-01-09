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
comp.wy <- comp.wy %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
cwm.wy <- read.csv("data/cwm_wy.csv")# Wyoming CWM data
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities
# combine to master df (remove spp columns for now)
allwy <- merge(comp.wy[,-c(11:66)],cwm.wy, by=c("year","trt","block","subplot"))
allwy$trt <- as.factor(allwy$trt)

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

#vizualize
# WY
ggplot(allwy, aes(y=nativecov, x=leafn, lty=drought.y, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
allwy$trt <- relevel(allwy$trt, ref = "rand") #make random communities the reference level
anova(lm(nativecov~leafn*drought.y*trt*year, data=allwy))
ggplot(allwy, aes(y=plot.tveg, x=srl, shape=drought.y, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
anova(lm(nativecov~srl*drought.y*trt*year, data=allwy))
library(lme4)
anova(lmer(nativecov~srl*drought.y*trt+(1|year), data=allwy))
library(lmerTest)
summary(lmer(sqrt(nativecov)~srl+drought.y+trt+year+(1|plot), data=allwy))
anova(lmer(sqrt(nativecov)~srl+drought.y+trt+year+(1|plot), data=allwy))
summary(lmer(sqrt(nativecov)~srl*drought.y*trt*year+(1|plot), data=allwy))

ggplot(allwy, aes(y=plot.tveg, x=ldmc, shape=drought.y, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
ggplot(allwy, aes(y=plot.tveg, x=lop, shape=drought.y, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
hist(allwy$ldmc)

ggplot(allwy, aes(y=plot.tveg, x=rootdiam, shape=drought.y, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
ggplot(allwy, aes(y=plot.tveg, x=veg, shape=drought.y, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)


ggplot(allwy, aes(y=nativecov, x=leafn, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought.y~year)
ggplot(allwy, aes(y=nativecov, x=srl, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought.y~year)




##testing dissim
# Determine how far off from the targets we were for EVERY treatment (refined in next steps)
allwy$dist_leafn <- abs(allwy$leafn - quantile(traits.wy$leafn,.25))
allwy$dist_srl <- abs(allwy$srl - quantile(traits.wy$srl,.7557))
allwy$dist_ldmc <- abs(allwy$ldmc - quantile(traits.wy$ldmc,.75))
allwy$dist_lop <- abs(allwy$lop - quantile(traits.wy$lop,.25))

# Calculate how far we were from the dt, ir, and fd objective. fd objective is simply highest diversity after normalization.
allwy$dist_dt <- rowSums(cbind(allwy$dist_ldmc,allwy$dist_lop))
allwy$dist_ir <- rowSums(cbind(allwy$dist_leafn,allwy$dist_srl))
allwy_roaq <- read.csv("data/cwm_raoq_wy.csv")
allwy$raoq <- allwy_roaq$full

ggplot(allwy, aes(y=sqrt(nativecov), x=dist_ldmc, col = trt))+ #plot
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(drought.y~year)
