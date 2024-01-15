#### CWM analysis, CA 
#### how do CWM's correlate with native cover (repeat drought models including CWM and CWM_dist and compare), 
#### with invasive cover (repeat invasion models including CWM and CWM_dist and compare), 
#### with FD, and ecosystem services 
## CWM as a predictor of native cover, invasive cover, and ecosystem services

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## load in data, clean and modify columns
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
cwm.ca$year <- as.factor(cwm.ca$year)
#make new sequence column
cwm.ca <- cwm.ca %>% mutate(yrorder = ifelse(year=="2021","1",
                                      ifelse(year=="2022","2",
                                      ifelse(year=="2023","3","0"))))
cwm.ca$yrorder <- as.numeric(cwm.ca$yrorder)
#add plot ID column (but give NA to target/predicted communities)
cwm.ca <- cwm.ca %>% 
  mutate(plot = paste(block, trt, sep = "."))
cwm.ca <- cwm.ca %>% filter(year != "0") #remove predicted/target communities
# combine to master df (remove spp columns for now)
allca <- merge(comp.ca[,-c(13:66)],cwm.ca, by.y=c("year","trt","block"), by.x=c("Year","trt","block"))
allca$trt <- as.factor(allca$trt)

# also combine CWM_distances dataframe to master df 
cadist <- read.csv("data/cwm_distances_ca.csv")
cadist <- cadist %>% select(-X) #%>% filter(trt.b.y!="target")
#break apart distances ID to make wider and merge together
cadist <- separate(cadist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist <- cadist %>% filter(trt!="target")
#cadist <- cadist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allca <- merge(allca,cadist, by.y=c("year","trt","block"), by.x=c("Year","trt","block"), all=T)
allca$trt <- as.factor(allca$trt)

# combine to master df (remove spp columns for now)
#allca <- merge(comp.ca[,-c(18:56)],cwm.ca, by=c("year","trt","block","water"))


#### Were CWM of DT communities correlated
#### I NEED TO calculate CWM at plot level for this model (?)
m1.wy.plot <- lmer(sqrt(plot.tveg) ~ trt * drought.x + (1 | year) + (1|plot), data = allwy)
allwy2 <- allwy %>% filter(subplot=="n") #only one subplot 
m2.wy.plot <- lmer(sqrt(plot.tveg) ~ drought.x * ldmc * lop * rootdiam + (1 | year) + (1|plot), data = allwy)
anova(m2.wy.plot)

allwy$trt <- relevel(allwy$trt, ref = "rand") #make random communities the reference level

m1.ca <- lmer(sqrt(native.cover) ~ trt * water.x + (1 | year) + (1|block), data = allca)
summary(m1.ca)
anova(m1.ca)

m2ldmc.wy <- lmer(sqrt(nativecov) ~ drought.x * LMA + (1 | year) + (1|block), data = allwy)
summary(m2ldmc.wy)
anova(m2ldmc.wy)
anova(m1.wy,m2ldmc.wy) #better than seeding trt
m2lop.wy <- lmer(sqrt(nativecov) ~ drought.x * seed.mass + (1 | year) + (1|block), data = allwy)
summary(m2lop.wy)
anova(m2lop.wy)
anova(m1.wy,m2lop.wy) #better than seeding trt
anova(m2ldmc.wy,m2lop.wy) #not better than ldmc
m2rd.wy <- lmer(sqrt(nativecov) ~ drought.x * rootdiam + (1 | year) + (1|block), data = allwy)
summary(m2rd.wy)
anova(m2rd.wy)
anova(m1.wy,m2rd.wy) #better than seeding trt
anova(m2ldmc.wy,m2rd.wy) #not better than ldmc
anova(m2lop.wy,m2rd.wy) #not better than lop
### Multivariate traits model/ value?

m3.wy <- lmer(sqrt(nativecov) ~ water.x * trt * dist + (1 | year) + (1|block), data = allwy)
summary(m3.wy)
anova(m3.wy)
anova(m1.wy,m3.wy) #better than seeding trt alone
anova(m2ldmc.wy,m3.wy) #better than CWM's alone
#when not including trt w/ dist, dist does not effect cover, but is interacting with drought


#### Were CWM of IR communities correlated
allwy23 <- allwy %>% filter(year=="2023")
allwy23$trt <- relevel(allwy23$trt, ref = "rand") # random as reference community

testwy <- allwy %>% mutate(log.brte = log(BRTE)) %>% 
  mutate(log.brte = ifelse(log.brte == -Inf, NA, log.brte))
m1.wy <- lmer(log.brte ~ trt * nativecov + (1 | block), data = allwy) #not working? random effects?
summary(m1.wy)
m2.wy <- lmer(log.brte ~ trt * nativecov * drought + (1 | block), data = allwy) #not working? random effects?
summary(m2.wy)
anova(m2.wy)

m1.wy.plot <- lmer(sqrt(plot.tveg) ~ trt * drought.x + (1 | year) + (1|plot), data = allwy)
allwy2 <- allwy %>% filter(subplot=="n") #only one subplot 
m2.wy.plot <- lmer(sqrt(plot.tveg) ~ drought.x * ldmc * lop * rootdiam + (1 | year) + (1|plot), data = allwy)
anova(m2.wy.plot)

allwy$trt <- relevel(allwy$trt, ref = "rand") #make random communities the reference level

m1.wy <- lmer(sqrt(nativecov) ~ trt * drought.x + (1 | year) + (1|block), data = allwy)
summary(m1.wy)
anova(m1.wy)
# m2.wy <- lmer(sqrt(nativecov) ~ drought.x * ldmc * lop * rootdiam + (1 | year) + (1|block), data = allwy)
# summary(m2.wy)
# anova(m2.wy)
m2ldmc.wy <- lmer(sqrt(nativecov) ~ drought.x * ldmc + (1 | year) + (1|block), data = allwy)
summary(m2ldmc.wy)
anova(m2ldmc.wy)
anova(m1.wy,m2ldmc.wy) #better than seeding trt
m2lop.wy <- lmer(sqrt(nativecov) ~ drought.x * lop + (1 | year) + (1|block), data = allwy)
summary(m2lop.wy)
anova(m2lop.wy)
anova(m1.wy,m2lop.wy) #better than seeding trt
anova(m2ldmc.wy,m2lop.wy) #not better than ldmc
m2rd.wy <- lmer(sqrt(nativecov) ~ drought.x * rootdiam + (1 | year) + (1|block), data = allwy)
summary(m2rd.wy)
anova(m2rd.wy)
anova(m1.wy,m2rd.wy) #better than seeding trt
anova(m2ldmc.wy,m2rd.wy) #not better than ldmc
anova(m2lop.wy,m2rd.wy) #not better than lop
### Multivariate traits model/ value?

m3.wy <- lmer(sqrt(nativecov) ~ drought.x * trt * dist + (1 | year) + (1|block), data = allwy)
summary(m3.wy)
anova(m3.wy)
anova(m1.wy,m3.wy) #better than seeding trt alone
anova(m2ldmc.wy,m3.wy) #better than CWM's alone
#when not including trt w/ dist, dist does not effect cover, but is interacting with drought




#### Figures
ggplot(allca, aes(y=native.cover, x=dist, col=water.x))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(Year~trt, scales="free")

ggplot(allca, aes(y=sqrt(dist), x=trt, fill=water.x))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  facet_wrap(~Year, scales="free")

allcasub <-allca %>% filter(Year=="2022")
ggplot(allcasub, aes(y=native.cover, x=dist, col=water.x))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(Year~trt, scales="free")
