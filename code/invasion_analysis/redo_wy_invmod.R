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
comp.wy <- read.csv("data/comp_wy_plot.csv") # Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
#correct factor columns
#comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy <- comp.wy %>% select(-drought)
#comp.wy$drought <- as.factor(comp.wy$drought)comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy <- comp.wy %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
# comp.wy$nativecov <- comp.wy$nativecov/100  # make native live veg % a proportion to match CA data
# comp.wy$totalcov <- comp.wy$totalcov/100  # make native live veg % a proportion to match CA data

cwm.wy <- read.csv("data/cwm_wy(plot).csv")# Wyoming CWM data
cwm.wy$year <- as.factor(cwm.wy$year)
#make new sequence column
cwm.wy <- cwm.wy %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.wy$yrorder <- as.numeric(cwm.wy$yrorder)
#add plot ID column (but give NA to target/predicted communities)
cwm.wy <- cwm.wy %>% 
  mutate(plot = paste(block, trt, sep = "."))
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities

cwm.wy <- cwm.wy %>% select(-c(rootdiam,veg)) #remove CWM rootdiam and veg column to avoid confusion
cwmFD <- read.csv("data/cwm_raoq_wy(plot).csv") #add FD for traits that need it (rootdiam/veg)
cwmFD <- cwmFD %>% select(block,trt,year,drought,rootdiam,veg) #only columns we need
cwm.wy <- merge(cwm.wy,cwmFD, all.x=T)

# combine to master df (remove spp columns for now,but keep BRTE)
allwy <- merge(comp.wy[,-c(5:11,13:60)],cwm.wy, by=c("year","trt","block"))
allwy$trt <- as.factor(allwy$trt)

# also combine CWM_distances dataframe to master df 
wydist <- read.csv("data/cwm_maxdistances_WY(plot).csv")
wydist <- wydist %>% select(-X) #%>% filter(trt.b.sub.y!="target")
#break apart distances ID to make wider and merge together
wydist <- separate(wydist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

allwy <- merge(allwy,wydist, by=c("year","trt","block"), all.x=T)
allwy23 <- allwy %>% filter(year=="2023")
allwy23$trt <- as.factor(allwy23$trt)
allwy23$block <- as.factor(allwy23$block)
allwy23$year <- as.factor(allwy23$year)
#allwy23$subplot <- as.factor(allwy23$subplot)
allwy23$drought <- as.factor(allwy23$drought)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
subvalid <- allwy23 %>% group_by(block,trt,year) %>%
  mutate(propnative = nativecov.plot/totcov.plot*100) %>% 
  filter(propnative < 80)
table(subvalid$year)
181/512*100 # 35% total observation to remove for inv. models
allwy23 <- allwy23 %>% 
  mutate(propnative = nativecov.plot/totcov.plot*100)

## Ensure levels are correctly compared in models
allwy23$trt <- relevel(allwy23$trt, ref = "rand") #make random communities the reference level
allwy23$drought <- relevel(allwy23$drought, ref = "cntl") #make ambient precip the reference level
#for vizuals
droughtcols <- c("cntl"="skyblue", "drt"="tomato") #create variable for color
#make invcov column
allwy23 <- allwy23 %>% mutate(invcov = totcov.plot-nativecov.plot) #(already a proportion)
#allwy23$BRTE <- allwy23$BRTE/100  # make native live veg % a proportion to match CA data

### data summary
# look at response variable in each dataset
# WY
hist(allwy23$BRTE)
hist(log(allwy23$BRTE)) #better logged
allwy23 %>% filter(BRTE!="0") %>% n_distinct() #only 192/512 (-250 or so not seeded or found) have BRTE
allwy23 <- allwy23 %>% mutate(log.brte = log(BRTE)) %>% 
  mutate(log.brte = ifelse(log.brte == -Inf, NA, log.brte))
#check additional variables
hist(allwy23$invcov)
hist(log(allwy23$invcov)) #better logged
allwy23 %>% filter(invcov!="0") %>% n_distinct() #only 259/512 (-259 not seeded or found) have BTRE
allwy23 <- allwy23 %>% mutate(log.inv = log(invcov)) %>% 
  mutate(log.inv = ifelse(log.inv == -Inf, NA, log.inv))

#### remove all plots where CWM could not be validly calculated
suballwy <- allwy23 %>% filter(propnative >= 80)

suballwy <- suballwy %>% filter(invaded == "1")

suballwy <- suballwy %>% 
  mutate(row = case_when(
    block %in% 1:8 ~ 1,
    block %in% 9:16 ~ 2,
    block %in% 17:24 ~ 3,
    block %in% 25:32 ~ 4,
    block %in% 33:40 ~ 5,
    block %in% 41:48 ~ 6,
    block %in% 49:56 ~ 7,
    TRUE ~ 8
  ))

#### Linear models of native cover ~ treatments:
#### (unsebsetted data could be used for this model)
### Cover-only model (fullest model):
#m0.wy <- lmer(log.brte ~ trt * drought + nativecov + (1 | block) + (1|plot.x), data = suballwy)
#summary(m0.wy) # plot can account for plot effect
## Error, grouping levels do not have enough observations. 
## removing plot for now.

### current cover-only model:
m1.wy <- lm(log.inv ~ trt * drought * nativecov, data = suballwy)
summary(m1.wy) 
#anova(m0.wy,m1.wy) 
## now block does not improve model statistically, may want to retain to account for design,
## but removing for now since it encapsulates essentially 0 variance.
emm.ca <- emmeans(m1.ca, c("trt","drought"))
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
emm.ca <- emmeans(m2leafn.wy, c("leafn","drought"))
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
anova(m1.wy,m2leafn.wy) #NOT better than model with seeding
#srl model
m2srl.wy <- lm(log.brte ~ srl * drought + nativecov, data = suballwy)
summary(m2srl.wy) 
anova(m2srl.wy) #all sig.
#emm.ca <- emmeans(m2lop.wy, c("lop","drought"))
#pairs(emm.ca)
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
# emm.ca <- emmeans(m2veg.wy, c("rootdiam","drought"))
# pairs(emm.ca)
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
anova(m2srl.wy,m2veg.wy) #not better than srl alone?
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








#### Linear models of native cover ~ treatments:
#### (unsebsetted data could be used for this model)
### Cover-only model (fullest model):
#m0.wy <- lmer(log.brte ~ trt * drought + nativecov + (1 | block) + (1|plot.x), data = suballwy)
#summary(m0.wy) # plot can account for plot effect
## Error, grouping levels do not have enough observations. 
## removing plot for now.

### current cover-only model:
m1.wy <- lm(log.brte ~ trt * drought * nativecov, data = suballwy)
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
emm.ca <- emmeans(m2leafn.wy, c("leafn","drought"))
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
anova(m1.wy,m2leafn.wy) #NOT better than model with seeding
#srl model
m2srl.wy <- lm(log.brte ~ srl * drought + nativecov, data = suballwy)
summary(m2srl.wy) 
anova(m2srl.wy) #all sig.
#emm.ca <- emmeans(m2lop.wy, c("lop","drought"))
#pairs(emm.ca)
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
# emm.ca <- emmeans(m2veg.wy, c("rootdiam","drought"))
# pairs(emm.ca)
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
anova(m2srl.wy,m2veg.wy) #not better than srl alone?
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




#save some for presentation 1/19
#CWM models
#lma
leafn.p <- ggplot(suballwy, aes(y=nativecov,x=leafn,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="cover BROTEC")+
  theme_classic()+
  xlim(-.5,1.5) #+ #remove outlier?
#seedmass
srl.p <- ggplot(suballwy, aes(y=nativecov,x=srl,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="")+
  theme_classic()#+
facet_wrap(~Year)
#rootdiam
veg.p<- ggplot(suballwy, aes(y=nativecov,x=veg,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="")+
  theme_classic()#+
facet_wrap(~Year)

###for 1/19
x <- ggpredict(m1.wy,c("nativecov [all]","drought","trt"), type = "re") #all smooths lines
tiff("figures/drought models/invasion_mod_wy.tiff", res=400, height = 4,width =5, "in",compression = "lzw")
plot(x, alpha = .1, show_title = F)+
  #labs(y="cover BROTEC", x="cover native species", col="drought trt")+
  scale_color_manual(values=c("darkblue","red4"))+
  theme_classic()
dev.off()
capture.output(summary(m1.wy), file="test.doc") #cov effected by drought, and dist:drought int.
emm.ca <- emmeans(m3.ca, c("dist","trt","drought"))#,"dist"))
?pairs(emm.ca,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.ca,pairwise~"drought","trt")
contrast(emm.ca)
test(emmeans(m3.ca, pairwise ~drought*trt*dist, at = list(dist = 0)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = .36)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = 1)))

###CWM figs
tiff("figures/drought models/invasion_cwm_wy.tiff", res=400, height = 2,width =7, "in",compression = "lzw")
(leafn.p + srl.p + veg.p + plot_layout(guides = 'collect')) 
dev.off()
capture.output(anova(m2lma.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2sm.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2rd.ca)[,c(3,5,6)], file="test2.doc") #cov effected by drought, and dist:drought int.


## for 2/7
quantile(suballca$distir,.05) #find lower .05
quantile(suballca$distir,.95) #find upper .95
ggplot(suballca, aes(y=log.invg ,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  #facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="absolute cover invasive grass", fill="seed trt")+
  theme_ggeffects()
x <- ggpredict(m3.ca,c("trt","distir [c(.75,2.56)]","drought","native.cover")) 
plot(x, show.title = F)

x <- ggpredict(m3.ca,c("native.cover","distir [c(.75,2.56)]","trt","drought")) 
plot(x, show.title = F)


#rootdiam model
#suballca$rootdiam <- normalize(suballca$rootdiam) # normalize FD
distmod <- lmer(log.invg ~ distir + (1 | structure), data = suballca)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=inv.grass.cov,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  # facet_wrap(~year)+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
