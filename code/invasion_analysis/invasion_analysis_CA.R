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
comp.ca$year<-as.factor(comp.ca$year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
#comp.ca$trt <- tolower(comp.ca$trt) #make these lower to match cwm dataframe
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
               by=c("year","trt","block","water"))
allca$trt <- as.factor(allca$trt)

# also combine CWM_distances dataframe to master df 
cadist <- read.csv("data/cwm_maxdistances_ca.csv")
cadist <- cadist %>% select(-X) #%>% filter(trt.b.y!="target")
#break apart distances ID to make wider and merge together
cadist <- separate(cadist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist <- cadist %>% filter(trt!="target")

#cadist <- cadist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allca <- merge(allca,cadist, 
               by=c("year","trt","block"), all=T)
allca22 <- allca %>% filter(year=="2022")
allca22$trt <- as.factor(allca22$trt)
allca22$plot.y <- as.factor(allca22$plot.y)
allca22$year <- as.factor(allca22$year)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
test <- allca22 %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  filter(propnative < 80)
18/204*100 # only 9% total observation to remove for inv. models
allca22 <- allca22 %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100)

#make dought column
allca22 <- allca22 %>% mutate(drought = ifelse(water=="0.5","drt","cntl"))

## Ensure levels are correctly compared in models
allca22$trt <- relevel(allca22$trt, ref = "rand") # random as reference level
allca22$water <- relevel(allca22$water, ref = "1.25") # water as reference level (drought = treatment)
#for vizuals
droughtcols <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

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

# assess only invaded plots
allca22 <- allca22 %>% filter(fesper.seeded=="1")

#### remove all plots where CWM could not be validly calculated
suballca <- allca22 %>% filter(propnative >= 80)

#### Linear models of native cover ~ treatments:
### Cover-only model:
#m1.ca <- lmer(log.invg ~ trt * water * native.cover + (1 | plot.y), data = suballca)
## Error, grouping levels do not have enough observations. 
## removing plot for now since it encapsulates essentially 0 variance.

### current cover-only model:
### (unsubsetted data could be used for this model (no CWM))
m1.ca <- lmer(log.invg ~ trt * drought * native.cover + (1|structure), data = suballca)
summary(m1.ca) 
anova(m1.ca) #effect of native cover only
emm.ca <- emmeans(m1.ca, c("trt","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=native.cover,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
ggplot(suballca, aes(y=log.invg,x=trt,fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)#+
facet_wrap(~Year)
# library(visreg)
# visreg(m1.ca,"trt", by="water")
library(ggeffects)
x <- ggpredict(m1.ca,c("drought","trt","native.cover"),back_transform = T) #all smooths lines
plot(x, alpha = .1)
x <- ggpredict(m1.ca,c("native.cover","drought","trt"),back_transform = T) #all smooths lines
plot(x, alpha = .1)
### CWM trait model: (individual models only right now)
### These models show how CWM effected the relationship with an environmental
### variable across the whole site
#leaf N model
m2leafn.ca <- lmer(log.invg ~ drought * N * native.cover+ (1|structure), data = suballca)
summary(m2leafn.ca) 
anova(m2leafn.ca) # only water 
emm.ca <- emmeans(m2leafn.ca, c("N","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=N,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
# library(visreg)
# visreg(m2leafn.ca,"N", by="water")
#compare
anova(m1.ca,m2leafn.ca) #NOT better fit than model with seeding trt
#seed mass model
m2srl.ca <- lmer(log.invg ~ drought * SRL * native.cover+ (1|structure), data = suballca)
summary(m2srl.ca) 
anova(m2srl.ca) #SRL and native cover and water
emm.ca <- emmeans(m2srl.ca, c("SRL","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=SRL,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~trt)
x <- ggpredict(m2srl.ca,c("SRL","drought","native.cover"),back_transform = T) #all smooths lines
plot(x, alpha = .1)
# library(visreg)
# visreg(m2srl.ca,"SRL", by="water")
#compare
anova(m1.ca,m2srl.ca) #Not as good as seeding trt
anova(m2leafn.ca,m2srl.ca) #same as leafn?
#root mass fraction model
m2rmf.ca <- lmer(log.invg ~ drought * RMF * native.cover+ (1|structure), data = suballca)
summary(m2rmf.ca) 
anova(m2rmf.ca) #water and diam maybe kindof
emm.ca <- emmeans(m2rmf.ca, c("RMF","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=RMF,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
  facet_wrap(~trt)
# library(visreg)
# visreg(m2rmf.ca,"RMF", by="water")
#compare
anova(m1.ca,m2rmf.ca) #Not better fit than model with seeding
anova(m2leafn.ca,m2rmf.ca) #same as lma? not better?
anova(m2srl.ca,m2rmf.ca) #same as sm? not better?
### Multivariate traits model/ value?

### CWM distance model: 
## trt can be dropped to improve slightly
m3.ca <- lmer(log.invg ~ drought * distir * trt * native.cover + (1|structure), data = suballca)
summary(m3.ca) #bad model
anova(m3.ca) #only water
emm.ca <- emmeans(m3.ca, c("trt","drought","distir"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=log.invg,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)#+
facet_grid(Year~trt)
ggplot(suballca, aes(x=dist,y=native.cover,color=water))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
facet_wrap(~trt)
x <- ggpredict(m3.ca,c("distir","trt","native.cover","drought"),back_transform = T) #all smooths lines
plot(x, alpha = .1)
#compare
anova(m1.ca,m3.ca) #slightly better than seeding trt alone
#anova(m2leafn.ca,m3.ca) #all CWM model not better 
#water is highly correlated with out ability to hit our targets



#figures to easily view three-way interaction (one uses trt as x-axis, th other uses dist)
x <- ggpredict(m1.ca,c("dist [all]","trt","water")) #all smooths lines
ca_threeway <- plot(x, alpha = .1, show_title = F)+
  labs(y="cover native species", x="Euclidian distance from target", col="seeding trt")+
  theme_classic()
x <- ggpredict(m3.ca,c("trt","distir","drought")) #all smooths lines
plot(x)+
  theme_classic()
x <- ggpredict(m3.ca,c("trt","native.cover","drought")) #all smooths lines
plot(x)+
  theme_classic()


#save some for presentation 1/19
#CWM models
#lma
leafn.p <- ggplot(suballca, aes(y=log.invg,x=N,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="cover invasive grass")+
  theme_classic()+
  xlim(-.5,1.5) + #remove outlier?
#seedmass
srl.p <- ggplot(suballca, aes(y=log.invg,x=SRL,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="")+
  theme_classic()#+
facet_wrap(~Year)
#rootdiam
rmf.p<- ggplot(suballca, aes(y=log.invg,x=RMF,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="")+
  xlim(-.5,1)+
  theme_classic()#+
facet_wrap(~Year)

#alltogether
library(patchwork)
tiff("figures/drought models/drought_cwm_ca.tiff", res=400, height = 4,width =6, "in",compression = "lzw")
ca_threeway / (lma.p + sm.p + rd.p + plot_layout(guides = 'collect')) 
dev.off()


###for 1/19
library(ggeffects)
x <- ggpredict(m1.ca,c("native.cover [all]","drought"), type = "re") #all smooths lines
tiff("figures/drought models/invasion_mod_ca.tiff", res=400, height = 4,width =5, "in",compression = "lzw")
plot(x, alpha = .1, show_title = F)+
  labs(y="cover invasive grass", x="native species cover", col="drought trt")+
  scale_color_manual(values=c("darkblue","red4"))+
  theme_classic()
dev.off()
capture.output(anova(m1.ca)[,c(1,4,5)], file="test.doc") #cov effected by drought, and dist:drought int.
emm.ca <- emmeans(m1.ca, c("native.cover","trt","drought"))#,"dist"))
pairs(emm.ca,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.ca,pairwise~"drought","trt")
contrast(emm.ca)
test(emmeans(m3.ca, pairwise ~drought*trt*dist, at = list(dist = 0)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = .36)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = 1)))

###CWM figs
tiff("figures/drought models/invasion_cwm_ca.tiff", res=400, height = 2,width =7, "in",compression = "lzw")
(leafn.p + srl.p + rmf.p + plot_layout(guides = 'collect')) 
dev.off()
capture.output(anova(m2lma.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2sm.ca)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2rd.ca)[,c(3,5,6)], file="test2.doc") #cov effected by drought, and dist:drought int.



## for 2/7
anova(m3.ca)
quantile(suballca$distir,.05) #find lower .05
quantile(suballca$distir,.95) #find upper .95
summary(lm(log.invg~trt*drought, suballca))
inv.cov.ca.plot <- ggplot(suballca, aes(y=inv.grass.cov ,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  #facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="absolute cover invasive grass", fill="seed trt")+
  theme_ggeffects()
x <- ggpredict(m3.ca,c("trt","distir [c(.72,2.51)]","drought","native.cover")) 
plot(x, show.title = F)

x <- ggpredict(m3.ca,c("distir [c(.72,2.51)]","drought","trt", "native.cover")) 
 plot(x, show.title = F)

 library(emmeans)
emm.ca <- emmeans(m3.ca, c("native.cover","trt","drought","distir [c(.72,2.51)]"))#,"dist"))
 
 contrast(emm.ca, "consec")
 

#rootdiam model
#suballca$rootdiam <- normalize(suballca$rootdiam) # normalize FD
distmod <- lmer(log.invg ~ distir * drought + (1 | structure), data = suballca)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(suballca, aes(y=inv.grass.cov,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  stat_cor(label.y = 40)+ 
  stat_regline_equation(label.y = 45) +
 #  annotate("text", x = min(data$x), y = max(data$y), 
 #           label = paste("R^2 =", round(r_squared, 3), "p-value =", round(p_value, 3)),
 #           hjust = 0, vjust = 1)
 # #  stat_regline_equation(aes(label = paste(after_stat(eq.label), 
 #                      paste("R^2 =", format(summary(model)$r.squared, digits = 2)))))+#, 
 #                           # sep = "*\", \"p-value =", format(summary(model)$coefficients[8], digits = 3), sep = "*, \""))))+
  # facet_wrap(~year)+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

