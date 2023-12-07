#### cover analysis, are invasion resistent communities more resistent to invasion? 
#### (invasion by focal annual grasses (BRTE or FESPER) + volunteers)
#### 2023 is the only viable year for this analysis in WY. 
#### 2022 + 2023 had FESPER in CA but not reseeded in 2023(?) for 2022 only comparable year.

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## load in data, clean and modify columns
# WY
comp.wy <- read.csv("data/hpg_total.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year == "2023") #keep only 2023 data 
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
comp.wy <- comp.wy %>% unite(plot, c(block, trt, subplot), sep = ".", remove=F) # make unique plot variable
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
                                                         "nativecov","BG", "Litter","plot","invaded"), 
                                             names_from = "species", 
                                             values_from = "cover")
comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data
#could use early season BRTE measures, but not necessarily comparable with CA so not right now

# create log response ratio as response metric 
comp.wy.wide <- comp.wy.wide %>% group_by(block,trt) %>% mutate(diffnativecov = nativecov[invaded == 1] / nativecov[invaded == 0])


# CA
comp.ca <- read.csv("data/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)]
comp.ca <- comp.ca %>% filter(Year == "2022") #keep only 2022 data right now
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
comp.ca$structure <- as.factor(comp.ca$structure)

# create log response ratio as response metric 
test <- comp.ca %>% group_by(block,trt) %>% 
  filter(mono==0)%>%
  #filter(all(c(0, 1) %in% fesper.present)) %>%
  mutate(diffnativecov = comp.ca$native.cover[comp.ca$fesper.seeded == 1] / comp.ca$native.cover[comp.ca$fesper.seeded == 0])


## data summary
## look at response variable in each dataset
# WY
hist(comp.wy.wide$BRTE)
hist(log(comp.wy.wide$BRTE)) # better logged

# CA
hist(comp.ca$FESPER)
hist(log(comp.ca$FESPER)) #better logged
hist(comp.ca$inv.grass.cov)
hist(log(comp.ca$inv.grass.cov)) #better logged

hist(log(comp.ca$inv.grass.cov))
hist(log(test$diffnativecov))


### visualize
# WY
ggplot(comp.wy.wide, aes(y=log(BRTE), x=nativecov, col=drought))+ #adding drought trt?
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trt)
ggplot(comp.wy.wide, aes(y=log(BRTE), x=trt, fill=drought))+ #adding drought trt?
  geom_boxplot()+
  scale_fill_manual(values=c("blue","red"))+
  facet_wrap(~year)

ggplot(test, aes(y=log(diffnativecov), x=nativecov, col=drought))+ #adding drought trt?
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
ggplot(test, aes(y=log(diffnativecov), x=trt, fill=drought))+ #adding drought trt?
  geom_boxplot()+
  scale_fill_manual(values=c("blue","red"))+
  facet_wrap(~year)

# CA
ggplot(comp.ca, aes(y=log(FESPER), x=native.cover, col=water))+ #adding water trt?
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trt)
ggplot(comp.ca, aes(y=log(inv.grass.cov), x=native.cover, col=water))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~trt)
ggplot(comp.ca, aes(y=log(diffnativecov), x=trt, fill=drought))+ #adding drought trt?
  geom_boxplot()+
  scale_fill_manual(values=c("blue","red"))+
  facet_wrap(~year)

ggplot(comp.ca, aes(y=log(FESPER), x=trt, fill=water))+ #adding drought trt?
  geom_boxplot()+
  scale_fill_manual(values=c("red","blue"))+
  facet_wrap(~year)
ggplot(comp.ca, aes(y=log(inv.grass.cov), x=trt, fill=water))+ #adding drought trt?
  geom_boxplot()+
  scale_fill_manual(values=c("red","blue"))+
  facet_wrap(~year)

hist(log(comp.wy.wide$BRTE)) #fix BRTE to work here 

### modelling invasive cover as a function of seeding community + native cover OR drought trt
# WY 
comp.wy.wide$trt <- relevel(comp.wy.wide$trt, ref = "r") # random as reference community

testwy <- comp.wy.wide %>% mutate(log.brte = log(BRTE)) %>% 
  mutate(log.brte = ifelse(log.brte == -Inf, NA, log.brte))
m1.wy <- lmer(log.brte ~ trt * nativecov + (1 | block), data = test) #not working? random effects?
summary(m1.wy)
m2.wy <- lmer(log.brte ~ trt * nativecov * drought + (1 | block), data = test) #not working? random effects?
summary(m2.wy)
anova(m2.wy)


log(comp.wy.wide$BRTE)
#new with daniel
m1.wy <- lmer(BRTE ~ trt * nativecov + (1 | block), data = comp.wy.wide) #not working? random effects?
m1.wy <- lm(log(diffnativecov) ~ trt * drought * nativecov + (1|block), data = comp.wy.wide) 
summary(m1.wy)


# CA
comp.ca$trt <- relevel(comp.ca$trt, ref = "R")
comp.ca$water <- relevel(comp.ca$water, ref = "1.25")

testca <- comp.ca %>% mutate(log.fesper = log(FESPER)) %>% 
  mutate(log.fesper = ifelse(log.fesper == -Inf, NA, log.fesper))

m1.ca <- lm(log.fesper ~ trt * native.cover * water, data = testca) 
summary(m1.ca)
m2.ca <- lmer(log.fesper ~ trt * native.cover * water + (1|block), data = testca) 
summary(m2.ca)

testca <- comp.ca %>% mutate(log.invg = log(FESPER)) %>% 
  mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))

m1.ca <- lm(log.invg ~ trt * native.cover * water, data = testca) 
summary(m1.ca)
m2.ca <- lmer(log.invg ~ trt * native.cover * water + (1|block), data = testca) 
summary(m2.ca)


m1.ca <- lmer(FESPER ~ trt * native.cover + (1 | Year) + (1 | block), data = comp.ca)
m1.ca <- lm(FESPER ~ trt * native.cover, data = comp.ca) # BROMAD ande BRODIA even worse
summary(m1.ca)

