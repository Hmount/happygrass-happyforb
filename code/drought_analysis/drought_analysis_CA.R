#### Analysis of drought treatments, CA 
#### how does drought effect native species in the different seeded communities? with 
#### different community CWM traits? and/or with our distance from targets?
#### 2 Q's: Are DT communities more tolerant of drought than random? than FD? 

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
allca$trt <- as.factor(allca$trt)
allca$plot.y <- as.factor(allca$plot.y)
allca$Year <- as.factor(allca$Year)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
test <- allca %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  filter(propnative < 80)
table(test$Year)
table(comp.ca$year)
(9+18+57)/(210*3)*100 # only 13% total observation to remove
allca <- allca %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100)

#make dought column
allca <- allca %>% mutate(drought = ifelse(water=="0.5","drt","cntl"))

## Ensure levels are correctly compared in models
allca$trt <- relevel(allca$trt, ref = "rand") # random as reference level
allca$drought <- as.factor(allca$drought)
allca$drought <- relevel(allca$drought, ref = "cntl") # water as reference level (drought = treatment)

#for vizuals
droughtcols <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# CA
hist(allca$native.cover) # native cover in CA is not skewed
## I am transforming WY data for normality, doing so here to match, little
## to no differences in model from this change.
# hist(sqrt(allca$native.cover)) # transform to match, but does not need (all models match both ways)
allca <- allca %>% 
  mutate(nativecov_tran = sqrt(native.cover))

#### remove all plots where CWM could not be validly calculated
suballca <- allca %>% filter(propnative >= 80)

#### Linear models of native cover ~ treatments:
### Cover-only model:
#m0.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year) + (1 | plot.y), data = test)
## I am getting a singularity error (probably from very small variance in random effects)
## models runs with warning, but removing plot for now since it encapsulates 
## essentially 0 variance.

### current cover-only model: (extremely similar w/ sqrt)
### (unsebsetted data could be used for this model)
m1.ca <- lmer(sqrt(native.cover) ~ trt * drought + (1 | Year), data = suballca)
summary(m1.ca) 
anova(m1.ca) #effect of water
emm.ca <- emmeans(m1.ca, c("trt","drought"))
pairs(emm.ca,simple="trt")
#view
ggplot(suballca, aes(y=native.cover,x=trt,fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)#+
  facet_wrap(~Year)

### CWM trait model: (individual models only right now)
### These models show how CWM effected the relationship with an environmental
### variable across the whole site
#lma model
m2lma.ca <- lmer(sqrt(native.cover) ~ drought * LMA + (1 | Year), data = suballca)
summary(m2lma.ca) 
anova(m2lma.ca) #drought and lma important
emm.ca <- emmeans(m2lma.ca, c("drought","LMA"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=native.cover,x=LMA,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  xlim(-.8,1.5) #+ #remove outlier?
  facet_wrap(~Year)
#compare
anova(m1.ca,m2lma.ca) #Equivalent fit as model with seeding
#seed mass model
m2sm.ca <- lmer(sqrt(native.cover) ~ drought * seed.mass + (1 | Year), data = suballca)
summary(m2sm.ca) 
anova(m2sm.ca) #only drought
emm.ca <- emmeans(m2sm.ca, c("seed.mass","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=native.cover,x=seed.mass,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
  facet_wrap(~Year)
#compare
anova(m1.ca,m2sm.ca) #Not as good as seeding trt
anova(m2lma.ca,m2sm.ca) #same as lma? lma has lower AIC tho
#rootdiam model
suballca$rootdiam <- normalize(suballca$rootdiam) # normalize FD
m2rd.ca2 <- lmer(sqrt(native.cover) ~ drought * rootdiam + (1 | Year), data = suballca)
summary(m2rd.ca2) 
anova(m2rd.ca2) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=native.cover,x=rootdiam,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_wrap(~Year)
#compare
anova(m1.ca,m2rd.ca) #Equivalent fit as model with seeding
anova(m2lma.ca,m2rd.ca) #same as lma? lma has lower AIC tho
anova(m2sm.ca,m2rd.ca) #same as sm? sm has WAY lower AIC tho
### Multivariate traits model/ value?

### CWM distance model: 
hist(suballca$dist)
hist(sqrt(suballca$dist))
## trt can be dropped
m3.ca <- lmer(sqrt(native.cover) ~ trt * dist * drought + (1 | Year), data = suballca)
summary(m3.ca)
anova(m3.ca) #cov effected by drought, and dist:drought int.
emm.ca <- emmeans(m3.ca, c("dist","trt","drought"))#,"dist"))
pairs(emm.ca)#,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.ca,pairwise~"drought","trt")
contrast(emm.ca)

#view
ggplot(suballca, aes(y=native.cover,x=dist,color=drought))+
  #geom_point()+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcols)+
  facet_wrap(~trt,scales="free")+
  theme_classic()
#compare
anova(m1.ca,m3.ca) #Better than seeding trt alone!(more like the same as)
anova(m2lma.ca,m3.ca) #NOT better than/equivalent to LMA or rootdiam CWM
anova(m2sm.ca,m3.ca) #slightly better than CWM sm alone
#when not including trt w/ dist*drought is a strong model
#drought is highly correlated with out ability to hit our targets

#figures to easily view three-way interaction (one uses trt as x-axis, th other uses dist)
# #x <- ggpredict(m3.ca,c("dist [all]","drought","trt"), type="re"),"drought")),type = "re") #all smooths lines
# #ca_threeway <- plot(x, alpha = .1, show_title = F)+
#   labs(y="cover native species", x="Euclidian distance from target", col="seeding trt")+
#   theme_classic()
# x <- ggpredict(m3.ca,c("trt","dist","drought")) #all smooths lines
# plot(x)+
#   theme_classic()

ggplot(suballca, aes(y=native.cover,x=dist,color=trt))+
  #geom_point()+
  geom_smooth(method = "lm",alpha=.25)+
  #scale_color_manual(values=droughtcols)+
  facet_wrap(~drought)+
  theme_classic()
ggplot(suballca, aes(y=native.cover,x=dist,color=drought))+
  #geom_point()+
  geom_smooth(method = "lm",alpha=.25)+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)+
  theme_classic()
#save some for presentation 1/19
#CWM models
#lma
lma.p <- ggplot(suballca, aes(y=native.cover,x=LMA,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="cover native species")+
  theme_classic()+
  xlim(-.5,1.5) #+ #remove outlier?
#seedmass
sm.p <- ggplot(suballca, aes(y=native.cover,x=seed.mass,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="")+
  theme_classic()#+
facet_wrap(~Year)
#rootdiam
rd.p<- ggplot(suballca, aes(y=native.cover,x=rootdiam,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="")+
  theme_classic()#+
facet_wrap(~Year)

#alltogether
library(patchwork)
tiff("figures/drought models/drought_cwm_ca.tiff", res=400, height = 4,width =6, "in",compression = "lzw")
(lma.p + sm.p + rd.p + plot_layout(guides = 'collect')) 
dev.off()


###for 1/19
x <- ggpredict(m3.ca,c("dist [all]","drought","trt"), type = "re") #all smooths lines
tiff("figures/drought models/drought_mod_ca.tiff", res=400, height = 4,width =5, "in",compression = "lzw")
plot(x, alpha = .1, show_title = F)+
     labs(y="cover native species", x="Euclidian distance from target", col="drought trt")+
  scale_color_manual(values=c("darkblue","red4"))+
  theme_classic()
dev.off()
m3.ca <- lmer(sqrt(native.cover) ~ trt * dist * drought + (1 | Year), data = suballca)
capture.output(anova(m3.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
emm.ca <- emmeans(m3.ca, c("dist","trt","drought"))#,"dist"))
?pairs(emm.ca,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.ca,pairwise~"trt","drought")
contrast(emm.ca)
test(emmeans(m3.ca, pairwise ~drought*trt*dist, at = list(dist = 0)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = .36)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = 1)))

###CWM figs
tiff("figures/drought models/drought_cwm_ca.tiff", res=400, height = 2,width =7, "in",compression = "lzw")
(lma.p + sm.p + rd.p + plot_layout(guides = 'collect')) 
dev.off()
capture.output(anova(m2lma.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2sm.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2rd.ca)[,c(3,5,6)], file="test2.doc") #cov effected by drought, and dist:drought int.
