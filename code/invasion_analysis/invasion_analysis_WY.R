#### Analysis of invasion treatments, CA 
#### how are inv species growth effected by native species in the different seeded communities? 
#### with different community CWM traits? and/or with our distance from targets?
#### 2 Q's: Are IR communities more tolerant of drought than random? than FD? 
#### OR: Is the likelihood of being invaded effected by any of the treatments? 

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

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
comp.wy$totalcov <- comp.wy$totalcov/100  # make native live veg % a proportion to match CA data

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

cwm.wy <- cwm.wy %>% select(-c(rootdiam,veg)) #remove CWM rootdiam column to avoid confusion
cwmFD <- read.csv("data/cwm_raoq_wy.csv") #add FD for traits that need it (rootdiam/veg)
cwmFD <- cwmFD %>% select(block,trt,year,subplot,drought,rootdiam,veg) #only columns we need
cwm.wy <- merge(cwm.wy,cwmFD, all.x=T)

# combine to master df (remove spp columns for now)
allwy <- merge(comp.wy[,-c(5,12:16,18:67)],cwm.wy, by=c("year","trt","block","subplot"))
allwy$trt <- as.factor(allwy$trt)

# also combine CWM_distances dataframe to master df 
wydist <- read.csv("data/cwm_distances_WY.csv")
wydist <- wydist %>% select(-X) #%>% filter(trt.b.sub.y!="target")
#break apart distances ID to make wider and merge together
wydist <- separate(wydist, trt.b.sub.y, into = c("trt", "block", "subplot", "year"), sep = "\\.")
wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

allwy <- merge(allwy,wydist, by=c("year","trt","block","subplot"), all.x=T)
allwy23 <- allwy %>% filter(year=="2023")
allwy23$trt <- as.factor(allwy23$trt)
allwy23$block <- as.factor(allwy23$block)
allwy23$year <- as.factor(allwy23$year)
allwy23$subplot <- as.factor(allwy23$subplot)
allwy23$drought <- as.factor(allwy23$drought)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
test <- allwy23 %>% group_by(block,trt,year,subplot) %>%
  mutate(propnative = nativecov/totalcov*100) %>% 
  filter(propnative < 80)
table(test$year)
181/512*100 # 35% total observation to remove for inv. models
allwy23 <- allwy23 %>% 
  mutate(propnative = nativecov/totalcov*100)

## Ensure levels are correctly compared in models
allwy23$trt <- relevel(allwy23$trt, ref = "rand") #make random communities the reference level
allwy23$drought <- relevel(allwy23$drought, ref = "cntl") #make ambient precip the reference level
#for vizuals
droughtcols <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# WY
hist(allwy23$BRTE)
hist(log(allwy23$BRTE)) #better logged
allwy23 %>% filter(BRTE!="0") %>% n_distinct() #only 192/512 (-250 or so not seeded or found) have BRTE
#check additional variables
# # I am not convinced of this variable 
# hist(allwy23$invg)
# hist(log(allca22$inv.grass.cov)) #better logged
# allca22 %>% filter(inv.grass.cov!="0") %>% n_distinct() #only 87/204 (-100 not seeded or found) have FESPER
allwy23 <- allwy23 %>% mutate(log.brte = log(BRTE)) %>% 
  mutate(log.brte = ifelse(log.brte == -Inf, NA, log.brte))

#### remove all plots where CWM could not be validly calculated
suballwy <- allwy23 %>% filter(propnative >= 80)

#### Linear models of native cover ~ treatments:
#### (unsebsetted data could be used for this model)
### Cover-only model (fullest model):
m0.wy <- lmer(log.brte ~ trt * drought + nativecov + (1 | block) + (1|plot.x), data = suballwy)
summary(m0.wy) # plot can account for plot effect
## Error, grouping levels do not have enough observations. 
## removing plot for now.

### current cover-only model:
m1.wy <- lm(log.brte ~ trt * drought + nativecov, data = suballwy)
summary(m1.wy) 
#anova(m0.wy,m1.wy) 
## now block does not improve model statistically, may want to retain to account for design,
## but removing for now since it encapsulates essentially 0 variance.
emm.ca <- emmeans(m1.ca, c("trt","water"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=log.brte,x=nativecov,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_wrap(~trt)
ggplot(suballwy, aes(y=log.brte,x=trt,fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)#+
facet_wrap(~Year)

### CWM trait model: (individual models only right now)
### These models show how CWM effected the relationship with an environmental
### variable across the whole site
### current cover-only model:
#leafn model
m2leafn.wy <- lm(log.brte ~ leafn * drought * nativecov, data = suballwy)
summary(m2leafn.wy) 
anova(m2leafn.wy) #bad model
emm.ca <- emmeans(m2ldmc.wy, c("ldmc","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=log.brte,x=leafn,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_wrap(~year)
library(visreg)
visreg(m2leafn.wy,"leafn", by="drought")
#compare
anova(m1.wy,m2leafn.wy) #Equivalent fit as model with seeding
#srl model
m2srl.wy <- lm(log.brte ~ srl * drought + nativecov, data = suballwy)
summary(m2srl.wy) 
anova(m2srl.wy) #all sig.
emm.ca <- emmeans(m2lop.wy, c("lop","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=log.brte,x=srl,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_wrap(~Year)
library(visreg)
visreg(m2srl.wy,"srl", by="drought")
#compare
anova(m1.wy,m2srl.wy) #NOT better than seeding trt alone
anova(m2leafn.wy,m2srl.wy) #not better than leaf n
#veg spread model
suballwy$veg <- normalize(suballwy$veg) # normalize FD
m2veg.wy <- lm(log.brte ~ veg * drought * nativecov, data = suballwy)
summary(m2veg.wy) 
anova(m2veg.wy) #all sig.
emm.ca <- emmeans(m2veg.wy, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=log.brte,x=veg,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_wrap(~year)
library(visreg)
visreg(m2veg.wy,"veg", by="drought")
#compare
anova(m1.wy,m2veg.wy) #NOT better than seeding trt alone
anova(m2leafn.wy,m2veg.wy) #?
anova(m2srl.wy,m2veg.wy) #not better than srl alone
### Multivariate traits model/ value?


### CWM distance model: 
## trt can be dropped to improve slightly
m3.wy <- lm(log.brte ~ drought * dist * trt * nativecov, data = suballwy)
summary(m3.wy) #bad model
anova(m3.wy) #only water
emm.ca <- emmeans(m3.wy, c("trt","drought","dist"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=log.brte,x=dist,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_grid(Year~trt)
ggplot(suballwy, aes(x=dist,y=nativecov,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
#compare
anova(m1.wy,m3.wy) #NOT better than seeding trt alone, but close
anova(m2veg.wy,m3.wy) #Not better fit with any CWM 
