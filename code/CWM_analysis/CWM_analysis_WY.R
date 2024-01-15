#### CWM analysis, WY 
#### how do CWM's correlate with native cover (repeat drought models including CWM and CWM_dist and compare), 
#### with invasive cover (repeat invasion models including CWM and CWM_dist and compare), 
#### with FD, and ecosystem services 
## CWM as a predictor of native cover, invasive cover, and ecosystem services

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

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
comp.wy$nativecov <- comp.wy$nativecov/100  # make native live veg % a proportion to match CA data
comp.wy$plot.tveg <- comp.wy$plot.tveg/100  # make native live veg % a proportion to match CA data
cwm.wy <- read.csv("data/cwm_wy.csv")# Wyoming CWM data
cwm.wy$year <- as.factor(cwm.wy$year)
#make new sequence column
cwm.wy <- cwm.wy %>% mutate(yrorder = ifelse(year=="2021","1",
                                      ifelse(year=="2022","2",
                                      ifelse(year=="2023","3","0"))))
cwm.wy$yrorder <- as.numeric(cwm.wy$yrorder)
#add plot ID column (but give NA to target/predicted communities)
cwm.wy <- cwm.wy %>% 
  mutate(plot = ifelse(!is.na(subplot), paste(block, trt, subplot, sep = "."), NA))
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities
# combine to master df (remove spp columns for now)
allwy <- merge(comp.wy[,-c(5,11:66)],cwm.wy, by=c("year","trt","block","subplot"))
allwy$trt <- as.factor(allwy$trt)

# also combine CWM_distances dataframe to master df 
wydist <- read.csv("data/cwm_distances_WY.csv")
wydist <- wydist %>% select(-X) #%>% filter(trt.b.sub.y!="target")
#break apart distances ID to make wider and merge together
wydist <- separate(wydist, trt.b.sub.y, into = c("trt", "block", "subplot", "year"), sep = "\\.")
wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allwy <- merge(allwy,wydist, by=c("year","trt","block","subplot"), all.x=T)
allwy$trt <- as.factor(allwy$trt)


#HERE
#break apart distances ID to make wider and merge together
# wydist <- read.csv("data/cwm_distances_WY.csv")
# wydist <- separate(wydist, trt.b.sub.y, into = c("trt", "block", "subplot", "year"), sep = "\\.")
# wydist <- wydist %>% filter(trt!="target")
# #wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
# allwy <- merge(cwm.wy,wydist, by=c("year","trt","block","subplot","plot"), all.x=T)
# allwy$trt <- as.factor(allwy$trt)

allwy$drought <- as.factor(allwy$drought)

#set reference levels for modelling
allwy$trt <- relevel(allwy$trt, ref = "rand") #make random communities the reference level
allwy$drought <- relevel(allwy$drought, ref = "cntl") #make random communities the reference level


#### Were CWM of DT communities correlated
#### plot CWM are similar, models would match CA, but not a lot of reason to loose replication
#### when other methods are not aligned anyway
m1.wy <- lmer(sqrt(nativecov) ~ trt * drought + (1 | year) + (1|plot.x) + (1|block), data = allwy)
summary(m1.wy)
anova(m1.wy)

# m2.wy <- lmer(sqrt(nativecov) ~ drought.x * ldmc * lop * rootdiam + (1 | year) + (1|block), data = allwy)
# summary(m2.wy)
# anova(m2.wy)
m2ldmc.wy <- lmer(sqrt(nativecov) ~ drought * ldmc + (1 | year) + (1|plot.x) + (1|block), data = allwy)
summary(m2ldmc.wy)
anova(m2ldmc.wy)
anova(m1.wy,m2ldmc.wy) #NOT better than seeding trt
m2lop.wy <- lmer(sqrt(nativecov) ~ drought * lop + (1 | year) + (1|plot.x) + (1|block), data = allwy)
summary(m2lop.wy)
anova(m2lop.wy)
anova(m1.wy,m2lop.wy) #better than seeding trt
anova(m2ldmc.wy,m2lop.wy) #better than/same as ldmc
m2rd.wy <- lmer(sqrt(nativecov) ~ drought * rootdiam + (1 | year) + (1|plot.x) + (1|block), data = allwy)
summary(m2rd.wy)
anova(m2rd.wy)
anova(m1.wy,m2rd.wy) #better than seeding trt
anova(m2ldmc.wy,m2rd.wy) #better than/same as ldmc
anova(m2lop.wy,m2rd.wy) # better than/same as ldmc
### Multivariate traits model/ value?

m3.wy <- lmer(sqrt(nativecov) ~ drought * dist + (1 | year) + (1|plot.x) + (1|block), data = allwy)
summary(m3.wy)
anova(m3.wy)
anova(m1.wy,m3.wy) #WAY better than seeding trt alone
anova(m2rd.wy,m3.wy) #better than CWM's alone
#when not including trt w/ dist, remains having string effect and model is strong


#### Were CWM of IR communities correlated
allwy23 <- allwy %>% filter(year=="2023")
allwy23$trt <- relevel(allwy23$trt, ref = "rand") # random as reference community

allwy23 <- allca %>% filter(Year=="2022")
allwy23$trt <- relevel(allwy23$trt, ref = "rand") # random as reference community
allwy23$water <- relevel(allwy23$water, ref = "1.25") # water as reference level (drought = treatment)

hist(allca23$FESPER)
hist(log(allca23$FESPER)) #better logged
testca <- allca %>% mutate(log.invg = log(FESPER)) %>% 
  mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))
sub <- allca23 %>% filter(inv.grass.cov!="0")
hist(allca23$inv.grass.cov)
hist(log(allca23$inv.grass.cov)) #better logged
testca <- allca %>% mutate(log.invg = log(inv.grass.cov)) %>% 
  mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))

#m1.ca <- lmer(log.invg ~ trt * native.cover + (1 | Year), data = testca)
m1.ca <- lmer(log.invg ~ trt * native.cover * water + (1 | Year), data = testca)
summary(m1.ca) 
anova(m1.ca)
emm.ca <- emmeans(m1.ca, c("trt","water"))
pairs(emm.ca)


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
ggplot(allwy, aes(y=sqrt(nativecov), x=sqrt(dist), col=drought.x))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(year~trt, scales="free")
ggplot(allwy, aes(y=sqrt(dist), x=trt, fill=drought.x))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  facet_wrap(~year, scales="free")


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
