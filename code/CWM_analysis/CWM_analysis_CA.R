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
allca <- merge(comp.ca[,-c(18:53,55:66)],cwm.ca, by.y=c("year","trt","block","water"), by.x=c("Year","trt","block","water"))
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
allca$plot.y <- as.factor(allca$plot.y)
allca$Year <- as.factor(allca$Year)

allca$trt <- relevel(allca$trt, ref = "rand") # random as reference level
allca$water <- relevel(allca$water, ref = "1.25") # water as reference level (drought = treatment)

#for vizuals
droughtcols <- c("cntl"="skyblue", "drt"="red") #create variable for color

#### Were CWM of DT communities correlated
#m1.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year) + (1 | plot.y), data = allca)
m1.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year), data = allca)
summary(m1.ca) 
anova(m1.ca)
emm.ca <- emmeans(m1.ca, c("trt","water"))
pairs(emm.ca)
## I am getting a singularity error (probably from very small variance in random effects)
## models runs with warning, but consider removing plot for now since it encapsulates 
## essentially 0 variance

m2lma.ca <- lmer(sqrt(native.cover) ~ water * LMA + (1 | Year), data = allca)
summary(m2lma.ca) 
anova(m2lma.ca)
anova(m1.ca,m2lma.ca) #NOT better than seeding trt
m2sm.ca <- lmer(sqrt(native.cover) ~ water * seed.mass + (1 | Year), data = allca)
summary(m2sm.ca) 
anova(m2sm.ca)
anova(m1.ca,m2sm.ca) #NOT better than seeding trt
anova(m2lma.ca,m2sm.ca) #Better than lma (?)
m2rd.ca <- lmer(sqrt(native.cover) ~ trt * Rdiam + (1 | Year), data = allca)
summary(m2rd.ca) 
anova(m2rd.ca)
anova(m1.ca,m2rd.ca) #NOT better than seeding trt
anova(m2lma.ca,m2rd.ca) #Better than lma (?)
anova(m2sm.ca,m2rd.ca) #Better than seedmass (?)
### Multivariate traits model/ value?

m3.ca <- lmer(sqrt(native.cover) ~ trt * dist * water + (1 | Year), data = allca)
summary(m3.ca)
anova(m3.ca)
anova(m1.ca,m3.ca) #Better than seeding trt alone!
anova(m2lma.ca,m3.ca) #NOT better than LMA CWM
anova(m2lma.ca,m3.ca) #Better than CWM's alone (sm or rd)
#when not including trt w/ dist*water is a strong model
#water is highly correlated with out ability to hit our targets










#### Were CWM of IR communities correlated
allca22 <- allca %>% filter(Year=="2022")
allca22$trt <- relevel(allca22$trt, ref = "rand") # random as reference community
allca22$water <- relevel(allca22$water, ref = "1.25") # water as reference level (drought = treatment)

# look at FESPER variable and across site
hist(allca22$FESPER)
hist(log(allca22$FESPER)) #better logged (below if needed)
# allca22 <- allca22 %>% mutate(log.invg = log(FESPER)) %>% 
#   mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))
sub <- allca22 %>% filter(FESPER!="0") #only 47 of 204 have FESPER
table(allca22$fesper.present,allca22$fesper.seeded) #(row,column)
# everywhere FESPER was not seeded should be removed since it does not have a paired plot
# and will only inflate 0's. But it is present where it was not seeded (and so are other)
# non-natives and we want to capture that those communities were subseptable to invasion
## currently filtering for all communities FESPER was seeded (>0 = FESPER could invade and
## how much, 0 = plot resisted invasion pressure) AND communities that ended up invaded 
## (>0 = plot did not resist invasion, 0 = cannot assess resistence to invasion pressure 
## must remove).

## However, when considering that non-focal inv sp were left in 2022, all plots (theoretically)
## had the same pressure of background invasion (with some facing higher via seeding)
sub <- allca22 %>% filter(inv.grass.cov!="0")
sub <- allca22 %>% filter(fesper.present!="0")
sub <- allca22 %>% filter(fesper.seeded!="0")

sub <- allca22 %>% filter(fesper.seeded!="0"|inv.grass.cov!="0")

hist(sub$inv.grass.cov)
hist(log(sub$inv.grass.cov)) #better logged
sub <- sub %>% mutate(log.invg = log(inv.grass.cov)) %>% 
  mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))

####TEST
mtest <- lmer(log.invg ~ trt * native.cover * water + (1|fesper.seeded), data = sub)
summary(mtest) 
anova(mtest)
emm.ca <- emmeans(mtest, c("trt","water"))
pairs(emm.ca)
#sig. less inv in ir.water than fd.drought
#sig. less inv in rand.water than fd.drought
#slightly sig. less inv in dt.water than fd.drought
#slightly sig. less inv in ir.water than ir.drought
### cover is not important and the model without it preforms better, consider removing
### unless to retain for matching CA/WY models. 

#m1.ca <- lmer(log.invg ~ trt * native.cover + (1 | Year), data = testca)
m1.ca <- lmer(log.invg ~ trt * native.cover * water + (1 | Year), data = testca)
summary(m1.ca) 
anova(m1.ca)
emm.ca <- emmeans(m1.ca, c("trt","water"))
pairs(emm.ca)
## model including models makes more since since drought treatment was also in effect and 
## definitely impacted communities

m2N.ca <- lmer(log.invg ~ N * native.cover * water + (1 | Year), data = testca)
summary(m2N.ca) 
anova(m2N.ca)
anova(m1.ca,m2N.ca) #NOT better than seeding trt
m2srl.ca <- lmer(log.invg ~ SRL * native.cover * water + (1 | Year), data = testca)
summary(m2srl.ca) 
anova(m2srl.ca)
anova(m1.ca,m2srl.ca) #NOT better than seeding trt
anova(m2N.ca,m2srl.ca) #Better than LMA (mostly the same?)
m2rmf.ca <- lmer(log.invg ~ RMF * native.cover * water + (1 | Year), data = testca)
summary(m2rmf.ca) 
anova(m2rmf.ca)
anova(m1.ca,m2rmf.ca) #NOT better than seeding trt
anova(m2N.ca,m2rmf.ca) #Better than lma (?)
anova(m2srl.ca,m2rmf.ca) #Better than seedmass (?)
### Multivariate traits model/ value?

m3.ca <- lmer(log.invg ~ trt * sqrt(dist) * native.cover * water + (1 | Year), data = testca)
summary(m3.ca)
anova(m3.ca)
anova(m1.ca,m3.ca) #Better than seeding trt alone
anova(m2N.ca,m3.ca) #NOT better than CWM's alone
anova(m2srl.ca,m3.ca) #slightly better than srl alone (?), but not if sqrt(dist)




#### Figures
ggplot(allca, aes(y=native.cover, x=dist, col=water))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(Year~trt, scales="free")

ggplot(allca, aes(y=sqrt(dist), x=trt, fill=water.x))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  facet_wrap(~Year, scales="free")

ggplot(testca, aes(y=native.cover, x=dist, col=water))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_grid(Year~trt, scales="free")

ggplot(allca22, aes(y=log.invg, x=native.cover, col=water))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trt, scales="free")
