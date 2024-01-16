#### Analysis of invasion treatments, CA 
#### how are inv species growth effected by native species in the different seeded communities? 
#### with different community CWM traits? and/or with our distance from targets?
#### 2 Q's: Are IR communities more tolerant of drought than random? than FD? 
#### OR: Is the likelihood of being invaded effected by any of the treatments? 

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

cwmFD <- read.csv("data/cwm_raoq_ca.csv") #add FD for traits that need it (rootdiam)
cwmFD <- cwmFD %>% select(block,trt,year,water,rootdiam) #only columns we need
cwm.ca <- merge(cwm.ca,cwmFD, all.x=T)
cwm.ca <- cwm.ca %>% select(-Rdiam) #remove CWM rootdiam column to avoid confusion

# combine to master df (remove spp columns for now)
allca <- merge(comp.ca[,-c(18:53,55:66)],cwm.ca,  #this merge drops monoculture plots
               by.y=c("year","trt","block","water"), 
               by.x=c("Year","trt","block","water"))
allca$trt <- as.factor(allca$trt)

# also combine CWM_distances dataframe to master df 
cadist <- read.csv("data/cwm_distances_ca.csv")
cadist <- cadist %>% select(-X) #%>% filter(trt.b.y!="target")
#break apart distances ID to make wider and merge together
cadist <- separate(cadist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist <- cadist %>% filter(trt!="target")

#cadist <- cadist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allca <- merge(allca,cadist, 
               by.y=c("year","trt","block"), 
               by.x=c("Year","trt","block"), all=T)
allca22 <- allca %>% filter(Year=="2022")
allca22$trt <- as.factor(allca22$trt)
allca22$plot.y <- as.factor(allca22$plot.y)
allca22$Year <- as.factor(allca22$Year)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
test <- allca22 %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  filter(propnative < 80)
18/204*100 # only 9% total observation to remove for inv. models
allca22 <- allca22 %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100)

## Ensure levels are correctly compared in models
allca22$trt <- relevel(allca22$trt, ref = "rand") # random as reference level
allca22$water <- relevel(allca22$water, ref = "1.25") # water as reference level (drought = treatment)
#for vizuals
droughtcols <- c("1.25"="skyblue", "0.5"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# CA
hist(allca22$FESPER)
hist(log(allca22$FESPER)) #better logged
allca22 %>% filter(FESPER!="0") %>% n_distinct() #only 47/204 (-100 not seeded or found) have FESPER
# I am not convinced of this variable 
hist(allca22$inv.grass.cov)
hist(log(allca22$inv.grass.cov)) #better logged
allca22 %>% filter(inv.grass.cov!="0") %>% n_distinct() #only 87/204 (-100 not seeded or found) have FESPER
allca22 <- allca22 %>% mutate(log.invg = log(inv.grass.cov)) %>% 
  mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))

#### remove all plots where CWM could not be validly calculated
suballca <- allca22 %>% filter(propnative >= 80)

#### Linear models of native cover ~ treatments:
### Cover-only model:
#m1.ca <- lmer(log.invg ~ trt * water * native.cover + (1 | plot.y), data = suballca)
## Error, grouping levels do not have enough observations. 
## removing plot for now since it encapsulates essentially 0 variance.

### current cover-only model:
### (unsubsetted data could be used for this model (no CWM))
m1.ca <- lm(log.invg ~ trt * water * native.cover, data = suballca)
summary(m1.ca) 
anova(m1.ca) #effect of water only
emm.ca <- emmeans(m1.ca, c("trt","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=native.cover,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
  facet_wrap(~trt)
ggplot(suballca, aes(y=log.invg,x=trt,fill=water))+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)#+
facet_wrap(~Year)

### CWM trait model: (individual models only right now)
### These models show how CWM effected the relationship with an environmental
### variable across the whole site
#leaf N model
m2leafn.ca <- lm(log.invg ~ water * N * native.cover, data = suballca)
summary(m2leafn.ca) 
anova(m2leafn.ca) # only water important
emm.ca <- emmeans(m2leafn.ca, c("N","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=N,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
library(visreg)
visreg(m2leafn.ca,"N", by="water")
#compare
anova(m1.ca,m2leafn.ca) #NOT better fit than model with seeding trt
#seed mass model
m2srl.ca <- lm(log.invg ~ water * SRL * native.cover, data = suballca)
summary(m2srl.ca) 
anova(m2srl.ca) #SRL and native cover
emm.ca <- emmeans(m2srl.ca, c("SRL","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=SRL,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
library(visreg)
visreg(m2srl.ca,"SRL", by="water")
#compare
anova(m1.ca,m2srl.ca) #Not as good as seeding trt
anova(m2leafn.ca,m2srl.ca) #same as leafn?
#root mass fraction model
m2rmf.ca <- lm(log.invg ~ water * RMF * native.cover, data = suballca)
summary(m2rmf.ca) 
anova(m2rmf.ca) #water and diam maybe kindof
emm.ca <- emmeans(m2rmf.ca, c("RMF","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=RMF,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
  facet_wrap(~trt)
library(visreg)
visreg(m2rmf.ca,"RMF", by="water")
#compare
anova(m1.ca,m2rmf.ca) #Not better fit than model with seeding
anova(m2leafn.ca,m2rmf.ca) #same as lma? not better?
anova(m2srl.ca,m2rmf.ca) #same as sm? not better?
### Multivariate traits model/ value?

### CWM distance model: 
## trt can be dropped to improve slightly
m3.ca <- lm(log.invg ~ water * dist * trt * native.cover, data = suballca)
summary(m3.ca) #bad model
anova(m3.ca) #only water
emm.ca <- emmeans(m3.ca, c("trt","water","dist"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=dist,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_grid(Year~trt)
ggplot(suballca, aes(x=dist,y=native.cover,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
facet_wrap(~trt)
#compare
anova(m1.ca,m3.ca) #NOT better than seeding trt alone
anova(m2leafn.ca,m3.ca) #all CWM model not better (? p-value not printing?)
#water is highly correlated with out ability to hit our targets




#### Logistic models ####
# add column for all inv grass presence absence
suballca <- suballca %>% mutate(inv.present = as.factor(ifelse(inv.grass.cov > 0, "1","0")))

### current cover-only logistic model:
### (unsubsetted data could be used for this model (no CWM))
logm1.ca <- glm(inv.present ~ trt * water, data = suballca, family = "binomial")
summary(logm1.ca) 
anova(logm1.ca) #effect of water only
# emm.ca <- emmeans(m1.ca, c("trt","water"))
# pairs(emm.ca)
#for for plot
pred <- predict(logm1.ca, type="response")
pred_df <- data.frame(native.cover = suballca$native.cover, pred)
#view
ggplot(suballca, aes(y=inv.present,x=native.cover))+
  geom_smooth(data = pred_df, aes(x = native.cover, y = pred), lty = 1, method = "glm")+
  #binomial_smooth(method="glm")+
  #geom_smooth(data = pred_df, aes(x = native.cover, y = pred), colour = "blue", method = "binimoal") +
  scale_color_manual(values=droughtcols)#+
facet_wrap(~trt)
ggplot(suballca, aes(y=log.invg,x=trt,fill=water))+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)#+
facet_wrap(~Year)

### CWM trait model: (individual models only right now)
### These models show how CWM effected the relationship with an environmental
### variable across the whole site
#leaf N model
m2leafn.ca <- lm(log.invg ~ water * N * native.cover, data = suballca)
summary(m2leafn.ca) 
anova(m2leafn.ca) # only water important
emm.ca <- emmeans(m2leafn.ca, c("N","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=N,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
library(visreg)
visreg(m2leafn.ca,"N", by="water")
#compare
anova(m1.ca,m2leafn.ca) #NOT better fit than model with seeding trt
#seed mass model
m2srl.ca <- lm(log.invg ~ water * SRL * native.cover, data = suballca)
summary(m2srl.ca) 
anova(m2srl.ca) #SRL and native cover
emm.ca <- emmeans(m2srl.ca, c("SRL","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=SRL,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
library(visreg)
visreg(m2srl.ca,"SRL", by="water")
#compare
anova(m1.ca,m2srl.ca) #Not as good as seeding trt
anova(m2leafn.ca,m2srl.ca) #same as leafn?
#root mass fraction model
m2rmf.ca <- lm(log.invg ~ water * RMF * native.cover, data = suballca)
summary(m2rmf.ca) 
anova(m2rmf.ca) #water and diam maybe kindof
emm.ca <- emmeans(m2rmf.ca, c("RMF","water"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=RMF,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_wrap(~trt)
library(visreg)
visreg(m2rmf.ca,"RMF", by="water")
#compare
anova(m1.ca,m2rmf.ca) #Not better fit than model with seeding
anova(m2leafn.ca,m2rmf.ca) #same as lma? not better?
anova(m2srl.ca,m2rmf.ca) #same as sm? not better?
### Multivariate traits model/ value?

### CWM distance model: 
## trt can be dropped to improve slightly
m3.ca <- lm(log.invg ~ water * dist * trt * native.cover, data = suballca)
summary(m3.ca) #bad model
anova(m3.ca) #only water
emm.ca <- emmeans(m3.ca, c("trt","water","dist"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=dist,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_grid(Year~trt)
ggplot(suballca, aes(x=dist,y=native.cover,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
#compare
anova(m1.ca,m3.ca) #NOT better than seeding trt alone
anova(m2leafn.ca,m3.ca) #all CWM model not better (? p-value not printing?)
#water is highly correlated with out ability to hit our targets
