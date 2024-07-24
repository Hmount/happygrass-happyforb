#### Analysis of drought treatments, WY 
#### how does drought effect native species in the different seeded communities? with 
#### different community CWM traits? and/or with our distance from targets?
#### 2 Q's: Are DT communities more tolerant of drought than random? than FD? 

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## load in data, clean and modify columns
# WY
comp.wy <- read.csv("data/comp_wy_plot.csv") # Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
#correct factor columns
#comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy <- comp.wy %>% select(-drought)
#comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
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
  mutate(plot = paste(block, trt, year, sep = "."))
# cwm.wy <- cwm.wy %>% 
#   mutate(plot = ifelse(!is.na(subplot), paste(block, trt, subplot, sep = "."), NA))
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities

cwm.wy <- cwm.wy %>% select(-c(rootdiam,veg)) #remove CWM rootdiam column to avoid confusion
cwmFD <- read.csv("data/cwm_raoq_wy(plot).csv") #add FD for traits that need it (rootdiam/veg)
cwmFD <- cwmFD %>% select(block,trt,year,drought,rootdiam,veg, full) #only columns we need
cwm.wy <- merge(cwm.wy,cwmFD, all.x=T)

# combine to master df (remove spp columns for now)
allwy <- merge(comp.wy[,-c(5:60)],cwm.wy, by=c("year","trt","block"))
allwy$trt <- as.factor(allwy$trt)

# also combine CWM_distances dataframe to master df 
wydist <- read.csv("data/cwm_maxdistances_wy(plot).csv")
wydist <- wydist %>% select(-X) #%>% filter(trt.b.sub.y!="target")
#break apart distances ID to make wider and merge together
wydist <- separate(wydist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

allwy <- merge(allwy,wydist, by=c("year","trt","block"), all.x=T)
allwy$trt <- as.factor(allwy$trt)
allwy$block <- as.factor(allwy$block)
allwy$year <- as.factor(allwy$year)
#allwy$subplot <- as.factor(allwy$subplot)
allwy$drought <- as.factor(allwy$drought)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
subvalid <- comp.wy %>% group_by(block,trt,year) %>%
  mutate(propnative = nativecov.plot/totcov.plot*100) %>% 
  filter(propnative < 80)
table(subvalid$year)
(17+132)/(512*3)*100 # only 10% total observation to remove
allwy <- allwy %>% 
  mutate(propnative = nativecov.plot/totcov.plot*100)

## Ensure levels are correctly compared in models
allwy$trt <- relevel(allwy$trt, ref = "rand") #make random communities the reference level
allwy$drought <- relevel(allwy$drought, ref = "cntl") #make ambient precip the reference level
#for vizuals
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# WY
hist(allwy$nativecov.plot) # native cover in WY is left skewed (extremely not normal)
hist(sqrt(allwy$nativecov.plot)) # good transformation for normalizing data
shapiro.test(sqrt(allwy$nativecov.plot)) #still not quite normal, but better
allwy <- allwy %>% 
  mutate(nativecov_tran = sqrt(nativecov.plot))

#### remove all plots where CWM could not be validly calculated
suballwy <- allwy %>% filter(propnative >= 80)

#### Linear models of native cover ~ treatments:
#### (unsebsetted data could be used for this model)
### Cover-only model (fullest model):
#m0.wy <- lmer(nativecov_tran ~ trt * drought + (1 | year) + (1 | block) + (1|plot), data = suballwy)
#summary(m0.wy) # plot can account for plot effect
### current cover-only model:
m1.wy <- lmer(nativecov_tran ~ trt * drought + (1 | year) + (1 | block), data = suballwy)
summary(m1.wy) 
#anova(m0.wy,m1.wy) 
## plot does not improve model statistically, may want to retain to account for design,
## but removing for now since it encapsulates very little variance.

### current cover-only model:
m1.wy <- lmer(nativecov_tran ~ trt * drought + (1 | year) + (1 | block), data = suballwy)
summary(m1.wy) # plot can account for plot effect
anova(m1.wy) #all effects sig.
emm.wy <- emmeans(m1.wy, c("trt","drought"))
pairs(emm.wy)
#view
ggplot(suballwy, aes(y=nativecov_tran,x=trt,fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)# +
  facet_wrap(~year)

### CWM trait model: (individual models only right now)
### These models show how CWM effected the relationship with an environmental
### variable across the whole site
#ldmc model
m2ldmc.wy <- lmer(nativecov_tran ~ ldmc * drought + (1 | year) + (1 | block), data = suballwy)
summary(m2ldmc.wy) 
anova(m2ldmc.wy) #all important important
emm.ca <- emmeans(m2ldmc.wy, c("ldmc","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=nativecov_tran,x=ldmc,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  xlim(-.8,1.5) + #remove outlier?
  facet_wrap(~year)
#compare
anova(m1.wy,m2ldmc.wy) #Equivalent fit as model with seeding
#lop model
m2lop.wy <- lmer(nativecov_tran ~ lop * drought + (1 | year) + (1 | block), data = suballwy)
summary(m2lop.wy) 
anova(m2lop.wy) #all sig.
emm.ca <- emmeans(m2lop.wy, c("lop","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=nativecov_tran,x=lop,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)
#compare
anova(m1.wy,m2lop.wy) #Better than seeding trt alone!
anova(m2ldmc.wy,m2lop.wy) #same as ldmc? lop has slightly lower AIC
#rootdiam model
suballwy$rootdiam <- normalize(suballwy$rootdiam) # normalize FD
m2rd.wy <- lmer(nativecov_tran ~ rootdiam * drought + (1 | year) + (1 | block), data = suballwy)
summary(m2rd.wy) 
anova(m2rd.wy) #all sig.
emm.ca <- emmeans(m2rd.wy, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=nativecov_tran,x=rootdiam,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)
#compare
anova(m1.wy,m2rd.wy) #Better than seeding trt alone!
anova(m2ldmc.wy,m2rd.wy) #same as ldmc? rd has slightly lower AIC
anova(m2lop.wy,m2rd.wy) #same as lop? rd has slightly lower AIC
### Multivariate traits model/ value?

suballwy$full <- normalize(suballwy$full) # normalize FD
m2rd.wy <- lmer(nativecov_tran ~ rootdiam * drought + (1 | year) + (1 | block), data = suballwy)
summary(m2rd.wy) 
anova(m2rd.wy) #all sig.
emm.ca <- emmea

### CWM distance model: 
## trt can be dropped
m3.wy <- lmer(nativecov_tran ~ trt * distdt * drought + (1 | year) + (1 | block), data = suballwy)
summary(m3.wy)
anova(m3.wy) #all sig! (except drought alone)
emm.wy <- emmeans(m3.wy, c("trt","drought","distdt") )
pairs(emm.wy,  at = list(distdt = 0.58))
contrast(emm.wy,simple=)

difflsmeans(m3.wy, at = list(distdt = 0.58))

#testing with daniel
glht(m3.wy, mcp(trt="Tukey", ))
multcomp::cdl(m3.wy, c("trt","drought","distdt"))
cld(m3.wy, c("trt","drought","distdt"))

#testing emmeasns
emm.wy <- emmeans(m3.wy, ~trt * distdt * drought)#,c("trt","drought","distdt"))
cld1 <- cld(emm.wy,adjust="sidak", Letters=letters)

### define coefficients of linear function directly
K <- diag(length(coef(m3.wy)))[-1,]
rownames(K) <- names(coef(m3.wy))[-1]
K 
### set up general linear hypothesis
glht(m3.wy, linfct = K)



#view
ggplot(suballwy, aes(y=nativecov,x=dist,color=trt))+
  #geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcolswy)+
facet_grid(~drought)
#compare
anova(m1.wy,m3.wy) #Better than seeding trt alone!
anova(m2rd.ca,m3.ca) #Better than and CWM trait alone! (except rootdiam is close)
#when not including trt w/ dist*water is a strong model
#water is highly correlated with out ability to hit our targets




#figures to easily view three-way interaction (one uses trt as x-axis, th other uses dist)
library(ggeffects)
x <- ggpredict(m3.wy,c("dist [all]","trt","drought")) #all smooths lines
wy_threeway <- plot(x, alpha = .1, show_title = F)+
  labs(y="cover native species", x="Euclidian distance from target", col="seeding trt")+
  theme_classic()
x <- ggpredict(m3.wy,c("trt","dist","drought")) #all smooths lines
plot(x)+
  theme_classic()

#save some for presentation 1/19
#CWM models
#lma
leafn.p <- ggplot(suballwy, aes(y=nativecov,x=leafn,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y="cover native species")+
  theme_classic()
  #xlim(-.5,1.5) #+ #remove outlier?
#seedmass
srl.p <- ggplot(suballwy, aes(y=nativecov,x=srl,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y="")+
  theme_classic()+
  facet_wrap(~year)
#rootdiam
veg.p<- ggplot(suballwy, aes(y=nativecov,x=rootdiam,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y="")+
  theme_classic()+
  facet_wrap(~year)

#alltogether
library(patchwork)
tiff("figures/drought models/drought_cwm_wy.tiff", res=400, height = 4,width =6, "in",compression = "lzw")
wy_threeway / (ldmc.p + lop.p + rd.p + plot_layout(guides = 'collect')) 
dev.off()

###for 1/19

#with dan and jen
m3.wy.y <- lmer(nativecov_tran ~ lop* ldmc*distdt * drought +(1|year)+ (1 | block), data = suballwy)

x <- ggpredict(m3.wy.y,c("ldmc [all]","drought","lop","distdt"), type = "re") #all smooths lines
x <- ggpredict(m3.wy,c("distdt [all]","drought","trt"), type = "re") #all smooths lines
plot(x, alpha = .1, show_title = F)+
  labs(y="cover native species", x="Euclidian distance from target", col="drought trt")+
  #scale_color_manual(values=c("darkblue","red4"))+
  #theme_classic()


x <- ggpredict(m3.wy,c("dist [all]","drought","trt"), type = "re") #all smooths lines
tiff("figures/drought models/drought_mod_wy.tiff", res=400, height = 4,width =5, "in",compression = "lzw")
plot(x, alpha = .1, show_title = F)+
  labs(y="cover native species", x="Euclidian distance from target", col="drought trt")+
  scale_color_manual(values=c("darkblue","red4"))+
  theme_classic()
dev.off()
m3.ca <- lmer(sqrt(native.cover) ~ trt * dist * drought + (1 | Year), data = suballca)
capture.output(anova(m3.wy)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
emm.ca <- emmeans(m3.wy, c("dist","trt","drought"))#,"dist"))
?pairs(emm.wy,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.wy,pairwise~"trt","drought")
contrast(emm.ca)
test<-(emmeans(m3.ca, pairwise ~drought*trt*distdt, at = list(distdt = c(0.58,1.58))))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = .36)))
test(emmeans(m3.ca, pairwise ~drought*dist, at = list(dist = 1)))

###CWM figs
tiff("figures/drought models/drought_cwm_wy.tiff", res=400, height = 2,width =7, "in",compression = "lzw")
(leafn.p + srl.p + veg.p + plot_layout(guides = 'collect')) 
dev.off()
capture.output(anova(m2ldmc.wy)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2lop.wy)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.
capture.output(anova(m2rd.wy)[,c(3,5,6)], file="test.doc") #cov effected by drought, and dist:drought int.


#### for 2/2 ###
tiff("figures/drought models/drought_mod_wy.tiff", res=400, height = 3,width =6, "in",compression = "lzw")
x <- ggpredict(m3.wy,c("trt","distdt","drought")) 
plot(x)+
  labs(x="seeding treatment", 
       y="absolute cover native species",
       title=" ")
#theme_classic()
dev.off()

#CWM's and FD
#CWM models
ldmc <- ggplot(suballwy, aes(y=nativecov.plot,x=ldmc,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y=" ")+
  theme(legend.position = "none")+
  #facet_wrap(~year)+
  theme_classic()
#xlim(-.5,1.5) #+ #remove outlier?
srl <- ggplot(suballwy, aes(y=nativecov.plot,x=lop, color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y="")+
  theme_classic()#+
  facet_wrap(~year)
rd <- ggplot(suballwy, aes(y=nativecov.plot,x=rootdiam,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y="", x="FD rootdiam")+
  theme_classic()#+
  facet_wrap(~year)
fd <- ggplot(suballwy, aes(y=nativecov.plot,x=full,color=drought))+
  geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  labs(y="", x="FD")+
  theme_classic()#+
  facet_wrap(~year)
#save
library(ggpubr)
tiff("figures/drought models/drought_CWM_wy.tiff", res=400, height = 4,width =8, "in",compression = "lzw")
annotate_figure(ggarrange(ldmc,srl,rd,fd, nrow=2, ncol=2, common.legend = T), 
                left ="absolute cover native species")
dev.off()

# for 2/5
#change quantiles
quantile(suballwy$distdt,.05) #find lower .05
quantile(suballwy$distdt,.95) #find upper .95

#make pairwise comparisons and letters to display
emm <- emmeans(m3.wy, specs = ~ trt + drought + distdt,  at = list(distdt = c(0.48,2.84)))
pairwise <- pairs(emm)
summary(pairwise, adjust = "tukey")
letters <- multcomp::cld(emm, Letters = LETTERS)


tiff("figures/drought models/drought_mod_wy.tiff", res=400, height = 3,width =6, "in",compression = "lzw")
x <- ggpredict(m3.wy,c("trt","distdt [c(.48,2.84)]","drought")) 

plot(x)+
  labs(x="seeding treatment", 
       y="absolute cover native species",
       title=" ")+
  scale_color_manual(values = c("1" = "slateblue", "2" = "lightsalmon3"), #orchid4 #olivdrab3
                     labels = c("1" = "0.48", "2" = "2.84"))#+
  #geom_text(data = letters, aes(x = trt, color=distdt, y = emmean, label = .group), hjust=2, vjust=-1.75, col="black") #add tukey labels

#theme_classic()
dev.off()

#test <- suballwy %>%
  mutate(distbin = ifelse(distdt <= 0.48,"close",
                          ifelse(distdt >= 2.84, "far", NA)))
#test <- test %>% filter(!is.na(distbin))
anova(lm(nativecov_tran~drought*trt*year, suballwy))
ggplot(suballwy, aes(y=nativecov_tran,x=trt))+
  geom_boxplot()+
  geom_line()+
  #geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  #scale_color_manual(values = c("close" = "slateblue", "far" = "lightsalmon3"))+
  facet_grid(drought~year,scales="free")#+
#theme_classic()


m3.wy <- lmer(nativecov_tran ~ trt*distdt * drought +(1|year)+ (1 | block), data = suballwy)
t <- lm(nativecov_tran~drought*trt*year, suballwy)
tuktest <- TukeyHSD(t)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought:year'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\2", rownames(letterstest)))
letterstest$year <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\3", rownames(letterstest)))
test <- allca %>% group_by(year, drought, trt) %>% summarise(yposition = quantile(dist,.8))
test <- merge(letterstest,test, by = c("year", "drought", "trt"))
test2 <- merge(test,allca, by = c("year", "drought", "trt"), all=T)

droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

ggplot(suballwy, aes(y=nativecov.plot,x=drought,col=trt))+
  geom_boxplot()+
  geom_line()+
  #geom_point(alpha=.3, pch=20)+
  geom_smooth(method = "lm")+
  # geom_text(aes(y=yposition,label = Letters), 
  #           position = position_dodge(width = 0.9), 
  #           vjust = -0.5,
  #           #angle = 15,
  #           size=3) +
  #scale_color_manual(values = c("close" = "slateblue", "far" = "lightsalmon3"))+
  facet_wrap(~year,scales="fixed")#+





#### 2/7/ comparison to loggr
ggplot(suballwy, aes(y=nativecov.plot,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="absolute cover native", fill="seed trt")+
  theme_ggeffects()
x <- ggpredict(m3.wy,c("trt","distdt [c(.56,1.56)]","drought", "year")) 
plot(x, show.title = F)

distmod <- lmer(nativecov_tran ~ distdt*trt + (1|year) + (1 | block), data = suballwy)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(suballwy, aes(y=nativecov.plot,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  labs(y="absolute cover native")+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

#together
ggarrange(gr.raw.fig.ca,gr.dist.fig.ca,gr.pred.fig.ca, nrow=3)


