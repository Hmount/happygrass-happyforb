#### Modelling Euclidean distances as a function of seeding treatment and drought 
#### to assess our ability to hit different targets in different communities and 
#### compare the success of those different targets to one another. 
#### We also examine the compositional dissimilarity (bray-curtis) as it relates 
#### to our trait-based targets in both sites.

# packages used
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)

#### read in data ####
cadist <- read.csv("data/cwm_maxdistances_ca.csv")# CA distance data
wydist <- read.csv("data/cwm_maxdistances_wy(plot).csv")# WY distance data

#### models and figures ####
####
#### CA
####

## seperate plot ID column into trt, block, and year
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
droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

## for each distance column create a model + figure showing distanmces to each target
## given each different seeding treatment
## Random 
summary(t <- aov(distr~trt*drought, cadist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- cadist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distr,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,cadist, by = c("drought", "trt"), all=T)
distr <- ggplot(test2, aes(y=distr, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolsca)+
  scale_alpha_manual(values = c(.4, .4, .4, 1))+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="exact CWM target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

## Functionally Diverse
summary(t <- aov(distfd~trt*drought, cadist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- cadist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distfd,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,cadist, by = c("drought", "trt"), all=T)
distfd <- ggplot(test2, aes(y=distfd, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolsca)+
  scale_alpha_manual(values = c(.4, 1, .4, .4))+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="FD target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

# Invasion resistant
summary(t <- aov(distir~trt*drought, cadist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- cadist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distir,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,cadist, by = c("drought", "trt"), all=T)
distir <- ggplot(test2, aes(y=distir, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolsca)+
  scale_alpha_manual(values = c(.4, .4, 1, .4))+
  #facet_wrap(~year, scales="fixed")+
  #coord_flip()+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="IR target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

#drought tolerant
summary(t <- aov(distdt~trt*drought, cadist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- cadist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distdt,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,cadist, by = c("drought", "trt"), all=T)
distdt <- ggplot(test2, aes(y=distdt, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolsca)+
  scale_alpha_manual(values = c(1, .4, .4, .4))+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="DT target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)
#attempt at adding yellow ellipse in R, but just do in ppt)
# annotation_custom(
#   grob = ellipseGrob(
#     x = unit(1, "inch"), 
#     y = unit(1, "inch"),
#     ar=1.5,
#     size=15,
#     angle=10,
#     gp = gpar(col = "yellow", fill = NA, lwd = 5)
#   ))#,
#   xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
# )
# annotate("point", x = 2, y = 3, size = 10, shape = 21, fill = "yellow", color = "yellow") +
# annotate("text", x = 2, y = 3.1, label = "Highlighted", color = "black")

## also create a model and figure comparing distances for each community to its intended
## target between seeding treatments
summary(t <- aov(targetdist~trt*drought, cadist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- cadist %>% group_by(drought, trt) %>% summarise(yposition = quantile(targetdist,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,cadist, by = c("drought", "trt"), all=T)
disttargetca <- ggplot(test2, aes(y=targetdist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolsca)+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="seeding treatment target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = c(.2,.9), legend.direction = "horizontal")+
  ylim(0,4)

library(ggpubr)
cadistplots <- ggarrange(distdt, distfd, distir, distr, ncol=2, nrow=2)
cadistanceplots<- ggarrange(cadistplots, disttargetca, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")



####
#### WY
####

## seperate plot ID column into trt, block, and year
wydist <- wydist %>% 
  separate(trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
wydist$block <- as.numeric(wydist$block)
wydist$year <- as.numeric(wydist$year)
## add drought column from cwm dataframe
cwm.wy <- read.csv("data/cwm_wy.csv")# California CWM data
wydrought <- cwm.wy %>% select(block,trt,year,drought) #drought and ID column
wydrought <- wydrought %>% distinct() #remove duplicates (from subplots)
wydist <- merge(wydist,wydrought, by=c("block", "trt","year")) #add drought treatment column
## set reference levels for modelling
wydist$trt <- as.factor(wydist$trt) #must be factor
wydist$trt <- relevel(wydist$trt, ref = "rand") #make random communities the reference level
wydist$drought <- as.factor(wydist$drought) #must be factor
wydist$drought <- relevel(wydist$drought, ref = "cntl") #make random communities the reference level
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

## for each distance column create a model + figure showing distanmces to each target
## given each different seeding treatment
## Random 
summary(t <- aov(distr~trt*drought, wydist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- wydist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distr,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,wydist, by = c("drought", "trt"), all=T)
distr <- ggplot(test2, aes(y=distr, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolswy)+
  scale_alpha_manual(values = c(.4, .4, .4, 1))+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="exact CWM target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

# Functionally diverse
summary(t <- aov(distfd~trt*drought, wydist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- wydist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distfd,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,wydist, by = c("drought", "trt"), all=T)
distfd <- ggplot(test2, aes(y=distfd, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolswy)+
  scale_alpha_manual(values = c(.4, 1, .4, .4))+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="FD target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

# Invasion resistant
summary(t <- aov(distir~trt*drought, wydist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- wydist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distir,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,wydist, by = c("drought", "trt"), all=T)
distir <- ggplot(test2, aes(y=distir, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolswy)+
  scale_alpha_manual(values = c(.4, .4, 1, .4))+
  #facet_wrap(~year, scales="fixed")+
  #coord_flip()+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="IR target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

#drought tolerant
summary(t <- aov(distdt~trt*drought, wydist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- wydist %>% group_by(drought, trt) %>% summarise(yposition = quantile(distdt,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,wydist, by = c("drought", "trt"), all=T)
distdt <- ggplot(test2, aes(y=distdt, x=trt, fill=drought, alpha=trt))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolswy)+
  scale_alpha_manual(values = c(1, .4, .4, .4))+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="DT target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = "none")+
  ylim(0,4)

## also create a model and figure comparing distances for each community to its intended
## target between seeding treatments
summary(t <- aov(targetdist~trt*drought, wydist))
tuktest <- TukeyHSD(t)
#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- wydist %>% group_by(drought, trt) %>% summarise(yposition = quantile(targetdist,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,wydist, by = c("drought", "trt"), all=T)
disttargetwy <- ggplot(test2, aes(y=targetdist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values = droughtcolswy)+
  #coord_flip()+
  #facet_wrap(~year, scales="fixed")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  labs(x=" ",y="seeding treatment target")+ #, fill="drought treatment")+
  theme_classic()+
  theme(legend.position = c(.2,.9), legend.direction = "horizontal")+
  ylim(0,4)

library(ggpubr)
wydistplots <- ggarrange(distdt, distfd, distir, distr, ncol=2, nrow=2)
wydistanceplot<- ggarrange(wydistplots, disttargetwy, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")



#### add/compare bray-Curtis dissimilarity to distance ####
## The relationship between our intended community composition and target CWM trait profile 
## should be tight and highly correlated.A poor relationship between these suggests that
## within a site communities differed from targets due to presence of unintended species in plots. 
## (as opposed to a fundamental issue in the premise of trait predictability/ intraspecific variation)
## See the "bray_curtis_calc.R script" for details on how these were calculated.

## read in bray-curtis data, create separated id columns, and make factors match 
bcdis.ca <- read.csv("data/bc_dissimilarity_ca.csv") #CA data
bcdis.ca <- bcdis.ca %>% 
  separate(id, into = c("trt", "block", "year"), sep = "\\.")
bcdis.ca$block <- as.numeric(bcdis.ca$block)
bcdis.ca$year <- as.numeric(bcdis.ca$year)
bcdis.ca <- bcdis.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cadist
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

## plot
#CA
bcplotca <- ggplot(cadisdist, aes(x=dist.ca, y=targetdist ,col=trt))+
  geom_point(pch=20)+
  geom_smooth(method="lm")+
  scale_color_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7) +
  #scale_color_viridis_d(option = "D", begin = 1, end = 0.1, alpha = 0.7) +
  labs(x=" ", 
       y="seeding treatment target",
       col="seeding trt")+
  theme_classic()
#WY
bcplotwy <- ggplot(wydisdist, aes(x=dist.wy, y=targetdist ,col=trt))+
  geom_point(pch=20)+
  geom_smooth(method="lm")+
  scale_color_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7) +
  #scale_color_viridis_d(option = "D", begin = 1, end = 0.1, alpha = 0.7) +
  labs(x=" ", 
       y="seeding treatment target",
       col="seeding trt")+
  theme_classic()

### possible arrangements, 

#option 1
## combine
# bcplots <- ggarrange(bcplotwy, bcplotca, ncol=2, nrow=1, common.legend = T)

#### combined figure for manuscript ####
# alldistanceplots <- ggarrange(wydistanceplots, cadistanceplots, bcplots, ncol=1, nrow=3)

wydistanceplot <- ggarrange(wydistplots, disttargetwy, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
cadistanceplot <- ggarrange(cadistplots, disttargetca, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
bcplots <- ggarrange(bcplotwy, bcplotca, ncol=2, nrow=1, common.legend = T)
p1 <-ggarrange(wydistanceplot, cadistanceplot,bcplots, nrow=3)
tiff("figures/test_figure.tiff", res=400, height = 8,width =6, "in",compression = "lzw")
p1
dev.off()

## option 2, best attempt yet, currently using this one (7/23/24):
wyp <- ggarrange(wydistplots, disttargetwy, ncol=2, nrow=1, common.legend = T)
cap <- ggarrange(cadistplots, disttargetca, ncol=2, nrow=1, common.legend = T)
p1 <- ggarrange(wyp, cap, ncol=1, nrow=2, common.legend = T)
p1 <- annotate_figure(p1,
                      bottom = text_grob("seeding treatment"))
bcplots <- ggarrange(bcplotwy, bcplotca, ncol=1, nrow=2, common.legend = T)
bcplots <- annotate_figure(bcplots, 
                      bottom = text_grob("Bray-Curtis dissimilarity"))
p2 <- ggarrange(p1, bcplots, ncol=2, nrow=1, widths = c(1,.5))
p2 <- annotate_figure(p2, 
                      left=text_grob("Euclidean distance from", rot=90))
p2
## export figure
tiff("figures/alldistances_figure.tiff", res=400, height = 7,width =9, "in",compression = "lzw")
p2
dev.off()
