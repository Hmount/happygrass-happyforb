#### Modelling Euclidean distances as a function of seeding treatment 
#### precipitation treatment to assess our ability to achieve specific 
#### trait-based targets for different restoration goals in communities. 
#### We also compare the success of communities with different seeding
#### treatments to one another. 
#### Last, we also examine the taxonomic dissimilarity (Bray-Curtis) as 
#### it relates to the distance between trait-based targets and realized
#### community traits in both sites.

#### updated 11/7/24 for simpler analysis + figures:
#### Initially, all communities were compared to each target (now commented 
#### out and below L412), but we are only interested in how each target 
#### compared to random controls. This reduced clutter in figure and better
#### highlights the important result of each figure.
#### specific changes included:
#### - random control is now RC
#### - removed/dropped the CWM in random comparison (this is not actually valid because
#### communities were chosen from log-normal distribution + not at all trait-based)

#### updated 1/23/25 to include year as a fixed effect in all models. Changes
#### are insignificant/marginal, but including year better accounts for the fact
#### we saw large interannual differences. 

## packages used
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)
library(ggpubr)
library(ggeffects)

## read in data
cadist <- read.csv("data/cwm_maxdistances_ca.csv")# CA distance data
wydist <- read.csv("data/cwm_maxdistances_wy(plot).csv")# WY distance data

## CA
## separate plot ID column into trt, block, and year
cadist <- cadist %>% 
  separate(trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist$block <- as.numeric(cadist$block)
cadist$year <- as.numeric(cadist$year)
## add drought column from cwm dataframe
cwm.ca <- read.csv("data/cwm_ca.csv")# California CWM data
cwm.ca <- cwm.ca %>% mutate(drought = ifelse(water=="0.5","drt","cntl")) #rename to dought column
cadrought <- cwm.ca %>% select(block,trt,year,drought) #drought and ID column
cadist <- merge(cadist,cadrought, by=c("block", "trt","year")) #add drought treatment column
## set reference levels for modelling
cadist$trt <- as.factor(cadist$trt) #must be factor
cadist$trt <- relevel(cadist$trt, ref = "rand") #make random communities the reference level
cadist$drought <- as.factor(cadist$drought) #must be factor
cadist$drought <- relevel(cadist$drought, ref = "cntl") #make random communities the reference level
cadist$year <- as.factor(cadist$year) #must be factor
droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

## WY
## seperate plot ID column into trt, block, and year
wydist <- wydist %>% 
  separate(trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
wydist$block <- as.numeric(wydist$block)
wydist$year <- as.numeric(wydist$year)
## add drought column from cwm dataframe
cwm.wy <- read.csv("data/cwm_wy.csv")# Wyoming CWM data
wydrought <- cwm.wy %>% select(block,trt,year,drought) #drought and ID column
wydrought <- wydrought %>% distinct() #remove duplicates (from subplots)
wydist <- merge(wydist,wydrought, by=c("block", "trt","year")) #add drought treatment column
## set reference levels for modelling
wydist$trt <- as.factor(wydist$trt) #must be factor
wydist$trt <- relevel(wydist$trt, ref = "rand") #make random communities the reference level
wydist$drought <- as.factor(wydist$drought) #must be factor
wydist$drought <- relevel(wydist$drought, ref = "cntl") #make random communities the reference level
wydist$year <- as.factor(wydist$year) #must be factor
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color


#### models and figures:

### CA
## Drought Tolerant
## are drought tolerant plots significantly closer to our DT target than random?
cadist.dt <- cadist %>% filter(trt=="dt"|trt=="rand")# subset for only DT and RC communities
summary(t <- aov(distdt~trt*drought+year, cadist.dt)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
dttemp <- cadist.dt %>% group_by(drought, trt) %>% summarise(yposition = quantile(distdt,.8))
dttemp <- merge(letters,dttemp, by = c("drought", "trt"))
dttemp <- merge(dttemp,cadist.dt, by = c("drought", "trt"), all=T)
#plot:
distdtca <- ggplot(dttemp, aes(y=distdt, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  guides(alpha = "none")+
  scale_x_discrete(labels = c("DT", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="... DT target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

## Invasion resistant 
## are invasion resistant plots significantly closer to our IR target than random?
cadist.ir <- cadist %>% filter(trt=="ir"|trt=="rand")# subset for only IR and RC communities
summary(t <- aov(distir~trt*drought+year, cadist.ir)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
irtemp <- cadist.ir %>% group_by(drought, trt) %>% summarise(yposition = quantile(distir,.8))
irtemp <- merge(letters,irtemp, by = c("drought", "trt"))
irtemp <- merge(irtemp,cadist.ir, by = c("drought", "trt"), all=T)
#plot:
distirca <- ggplot(irtemp, aes(y=distir, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  guides(alpha = "none")+
  scale_x_discrete(labels = c("IR", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="... IR target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

## Functionally Diverse
## are these plots significantly more functionally diverse (Rao) than random?
cadist.fd <- cadist %>% filter(trt=="fd"|trt=="rand")# subset for only FD and RC communities
summary(t <- aov(distfd~trt*drought+year, cadist.fd)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t) #run post-hoc for letters
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
fdtemp <- cadist.fd %>% group_by(drought, trt) %>% summarise(yposition = quantile(distfd,.8))
fdtemp <- merge(letters,fdtemp, by = c("drought", "trt"))
fdtemp <- merge(fdtemp,cadist.fd, by = c("drought", "trt"), all=T)
#plot:
distfdca <- ggplot(fdtemp, aes(y=distfd, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  guides(alpha = "none")+
  scale_x_discrete(labels = c("FD", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="... FD target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)


## Next, create a model and figure comparing distances for each community to 
## its intended target 
## Did seeding treatments differ in out ability to hit the target?
cadistsub <- cadist %>% filter(trt!="rand")
summary(t <- aov(targetdist~trt*drought+year, cadistsub))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
alltemp <- cadistsub %>% group_by(drought, trt) %>% summarise(yposition = quantile(targetdist,.8))
alltemp <- merge(letters,alltemp, by = c("drought", "trt"))
alltemp <- merge(alltemp,cadistsub, by = c("drought", "trt"), all=T)
disttargetca <- ggplot(alltemp, aes(y=targetdist, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca, labels=c("Ambient/ 
Addition", "Reduction"))+
  scale_x_discrete(labels = c("DT", "FD", "IR"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="Euclidean distance from seeding treatment target", fill="Precipitation 
treatment")+  theme_ggeffects()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  ylim(0,4)



### WY
## Drought Tolerant
## are drought tolerant plots significantly closer to our DT target than random?
wydist.dt <- wydist %>% filter(trt=="dt"|trt=="rand")# subset for only DT and RC communities
summary(t <- aov(distdt~trt*drought+year, wydist.dt)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
dttemp <- wydist.dt %>% group_by(drought, trt) %>% summarise(yposition = quantile(distdt,.8))
dttemp <- merge(letters,dttemp, by = c("drought", "trt"))
dttemp <- merge(dttemp,wydist.dt, by = c("drought", "trt"), all=T)
#plot:
distdtwy <- ggplot(dttemp, aes(y=distdt, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  scale_x_discrete(labels = c("DT", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="... DT target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

## Invasion resistant 
## are invasion resistant plots significantly closer to our IR target than random?
wydist.ir <- wydist %>% filter(trt=="ir"|trt=="rand")# subset for only IR and RC communities
summary(t <- aov(distir~trt*drought+year, wydist.ir)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
irtemp <- wydist.ir %>% group_by(drought, trt) %>% summarise(yposition = quantile(distir,.8))
irtemp <- merge(letters,irtemp, by = c("drought", "trt"))
irtemp <- merge(irtemp,wydist.ir, by = c("drought", "trt"), all=T)
#plot:
distirwy <- ggplot(irtemp, aes(y=distir, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  scale_x_discrete(labels = c("IR", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="... IR target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

## Functionally Diverse
## are these plots significantly more functionally diverse (Rao) than random?
wydist.fd <- wydist %>% filter(trt=="fd"|trt=="rand")# subset for only FD and RC communities
summary(t <- aov(distfd~trt*drought+year, wydist.fd)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t) #run post-hoc for letters
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
fdtemp <- wydist.fd %>% group_by(drought, trt) %>% summarise(yposition = quantile(distfd,.8))
fdtemp <- merge(letters,fdtemp, by = c("drought", "trt"))
fdtemp <- merge(fdtemp,wydist.fd, by = c("drought", "trt"), all=T)
#plot:
distfdwy <- ggplot(fdtemp, aes(y=distfd, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  scale_x_discrete(labels = c("FD", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.45), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="... FD target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)


## Next, create a model and figure comparing distances for each community to 
## its intended target 
## Did seeding treatments differ in out ability to hit the target?
wydistsub <- wydist %>% filter(trt!="rand")
summary(t <- aov(targetdist~trt*drought+year, wydistsub))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
letters$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letters)))
alltemp <- wydistsub %>% group_by(drought, trt) %>% summarise(yposition = quantile(targetdist,.8))
alltemp <- merge(letters,alltemp, by = c("drought", "trt"))
alltemp <- merge(alltemp,wydistsub, by = c("drought", "trt"), all=T)
disttargetwy <- ggplot(alltemp, aes(y=targetdist, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolsca)+
  scale_x_discrete(labels = c("DT", "FD", "IR"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            size=3) +
  labs(x=" ",y="Euclidean distance from seeding treatment target", fill="Drought treatment")+
  theme_ggeffects()+
  theme(legend.position = c(.3,.9), legend.direction = "vertical")+
  ylim(0,4)

### Only in Appendix S1 Figure S3 (see "pca.R" for construction)
### this relationship describes how well we did at making the experiment, not
### so much how well we did at making communities
## Last, plot bray-Curtis dissimilarity as it relates to Euclidean distance.
## The relationship between our intended community composition and target CWM trait profile
## should be tight and highly correlated. A tight relationship suggests that
## the replacement of species corresponds to the change in traits. But, a
## poor relationship between these suggests that the replacement of species did
## not lead to consistent changes in CWM, possibly suggesting high functional
## redundancy.
## See the "bray_curtis_calc.R script" for details on how these were calculated.

## read in bray-curtis data, create separated id columns, and make factors match
## CA
bcdis.ca <- read.csv("data/bc_dissimilarity_ca.csv") #CA data
bcdis.ca <- bcdis.ca %>%
  separate(id, into = c("trt", "block", "year"), sep = "\\.")
bcdis.ca$block <- as.numeric(bcdis.ca$block)
bcdis.ca$year <- as.numeric(bcdis.ca$year)
bcdis.ca <- bcdis.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cadist
## WY
bcdis.wy <- read.csv("data/bc_dissimilarity_wy.csv") #WY data
bcdis.wy <- bcdis.wy %>%
  separate(id, into = c("trt", "block", "year"), sep = "\\.")
bcdis.wy$block <- as.numeric(bcdis.wy$block)
bcdis.wy$year <- as.numeric(bcdis.wy$year)
bcdis.wy <- bcdis.wy %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cadist

## combine with Euclidean distance data (calculated above)
## (dis.site column = bray-curtis)
cadisdist <- merge(cadist, bcdis.ca)
wydisdist <- merge(wydist, bcdis.wy)

#model and plot:
## CA
summary(lm(targetdist~dist.ca*trt, cadisdist))
confint(lm(targetdist~dist.ca*trt, cadisdist))
bcplotca <- ggplot(cadisdist, aes(x=dist.ca, y=targetdist, col=trt))+
  geom_point(pch=20)+
  geom_smooth(method="lm")+
  scale_color_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7,
                        labels = c("dt" = "DT", "fd" = "FD",
                                   "ir" = "IR", "rand" = "RC")) +
  labs(x="Bray-Curtis dissimilarity",
       y="Euclidean distance to
seeding treatment target",
       col="Seeding
treatment")+
  theme_ggeffects()+
  theme(legend.position = "bottom")
#theme(legend.position = c(-.5,0), legend.direction = "horizontal")
## WY
summary(lm(targetdist~dist.wy*trt, wydisdist))
confint(lm(targetdist~dist.wy*trt, wydisdist))
bcplotwy <- ggplot(wydisdist, aes(x=dist.wy, y=targetdist ,col=trt))+
  geom_point(pch=20)+
  geom_smooth(method="lm")+
  scale_color_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7,
                        labels = c("dt" = "DT", "fd" = "FD",
                                   "ir" = "IR", "rand" = "RC"))+
  labs(x="Bray-Curtis dissimilarity",
       y="Euclidean distance to
seeding treatment target",
       col="Seeding
treatment")+
  theme_ggeffects()+
  theme(legend.position = "bottom")
#theme(legend.position = c(.3,.8))#, legend.direction = "horizontal")

## Combining
# get legends
drought_legend <- as_ggplot(get_legend(disttargetca))
#seed_legend <- as_ggplot(get_legend(bcplotca))

#first the distance plots
wydistplots <- ggarrange(distdtwy, distirwy, distfdwy, ncol=3, nrow=1, 
                         common.legend = T, legend = "none", 
                         labels = c("a","b","c"), hjust=c(-4,-4,-4))
# wydistbcplot <- ggarrange(wydistplots, bcplotwy, ncol=2, nrow=1,
#                           widths = c(1,.5), 
#                           legend = "none",
#                           labels = c(" ", "d"), hjust = -5)
wydistbcplot <- annotate_figure(wydistplots, right=text_grob("Wyoming", rot=270, hjust = .8),fig.lab.face="bold")

cadistplots <- ggarrange(distdtca, distirca, distfdca, ncol=3, nrow=1, 
                         common.legend = T, legend = "none", 
                         labels = c("d","e","f"), hjust=c(-4,-8,-4))
# cadistbcplot <- ggarrange(cadistplots, bcplotca, ncol=2, nrow=1,
#                           widths = c(1,.5), 
#                           legend = "none",
#                           labels = c(" ", "h"), hjust = -5)
cadistbcplot <- annotate_figure(cadistplots, right=text_grob("California", rot=270, hjust = .8), fig.lab.face="bold")

alldistbcplot <- ggarrange(wydistbcplot, cadistbcplot, ncol=1, nrow=2)
alldistbcplot <- annotate_figure(alldistbcplot, left=text_grob("Euclidean distance to...", rot=90))

legend <- ggarrange(drought_legend, nrow = 1)

alldistbcplot <- ggarrange(alldistbcplot, legend, nrow=2, heights = c(1.75,.25))
alldistbcplot

## export figure (add star and details in ppt.)
# tiff("figures/alldistances_figure.tiff", res=400, height = 7,width =9, "in",compression = "lzw")
# alldistbcplot
# dev.off()
tiff("figures/alldistances_figure.tiff", res=800, height = 5,width =7, "in") #,compression = "lzw")
alldistbcplot
dev.off()
