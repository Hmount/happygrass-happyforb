#### Analysis of invasion, WY 
#### how are inv species growth effected by native species in the different seeded communities? and
#### with our distance from targets?
#### Are IR communities or traits any more resistant to invasion than random? than FD? 
#### Are communities with more IR-target traits more tolerant of drought?
##(Is the likelihood of being invaded effected by any of the treatments?)

## packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggeffects)

## read in all data and clean to create master dataframe
comp.wy <- read.csv("data/comp_wy_plot.csv") # Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year == "2023") #remove pre-treatment data
#correct factor columns
comp.wy <- comp.wy %>% select(-drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy <- comp.wy %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

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
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities

cwm.wy <- cwm.wy %>% select(-c(rootdiam,veg)) #remove CWM rootdiam column to avoid confusion
cwmFD <- read.csv("data/cwm_raoq_wy(plot).csv") #add FD for traits that need it (rootdiam/veg)
cwmFD <- cwmFD %>% select(block,trt,year,drought,rootdiam,veg, full) #only columns we need
cwm.wy <- merge(cwm.wy,cwmFD, all.x=T)

# combine to master df (remove spp columns for now)
allwy23 <- merge(comp.wy[,-c(5:11,13:60)],cwm.wy, by=c("year","trt","block"))
allwy23$trt <- as.factor(allwy23$trt)

# also combine CWM_distances dataframe to master df 
wydist <- read.csv("data/cwm_maxdistances_wy(plot).csv")
#break apart distances ID to make wider and merge together
wydist <- separate(wydist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

allwy23 <- merge(allwy23,wydist, by=c("year","trt","block"), all.x=T)
allwy23$trt <- as.factor(allwy23$trt)
allwy23$block <- as.factor(allwy23$block)
allwy23$year <- as.factor(allwy23$year)
#allwy$subplot <- as.factor(allwy$subplot)
allwy23$drought <- as.factor(allwy23$drought)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
subvalid <- comp.wy %>% group_by(block,trt,year) %>%
  mutate(propnative = nativecov.plot/totcov.plot*100) %>% 
  filter(propnative < 80)
table(subvalid$year)
(17+132)/(512*3)*100 # only 10% total observation to remove
allwy23 <- allwy23 %>% 
  mutate(propnative = nativecov.plot/totcov.plot*100) %>% 
  mutate(propbrte = BRTE/totcov.plot) %>% #this one is not converted to percent
  mutate(propnativebrte = propnative+(BRTE*100))

## Ensure levels are correctly compared in models
allwy23$trt <- relevel(allwy23$trt, ref = "rand") #make random communities the reference level
allwy23$drought <- relevel(allwy23$drought, ref = "cntl") #make ambient precip the reference level
#for visuals
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# looking at plots as a whole even though we seeded into one of the subplots, so 256 total plots
# WY
hist(allwy23$BRTE)
hist(log(allwy23$BRTE)) #better logged
allwy23 %>% filter(BRTE!="0") %>% n_distinct() #only 192/512 (-250 or so not seeded or found) have BRTE (~75%)
allwy23 <- allwy23 %>% mutate(log.brte = log(BRTE)) %>% 
  mutate(log.brte = ifelse(log.brte == -Inf, NA, log.brte))
# #check additional variables
hist(allwy23$propbrte)
hist(log(allwy23$propbrte)) #better logged
allwy23 %>% filter(propbrte!="0") %>% n_distinct() #only 192/512 (-250 or so not seeded or found) have BRTE (~75%)
allwy23 <- allwy23 %>% mutate(log.propbrte = log(propbrte)) %>% 
  mutate(log.propbrte = ifelse(log.propbrte == -Inf, NA, log.propbrte))
# #make invcov column
# allwy23 <- allwy23 %>% mutate(invcov = totcov.plot-nativecov.plot) #(already a proportion)
# hist(allwy23$invcov)
# hist(log(allwy23$invcov)) #better logged
# allwy23 %>% filter(invcov!="0") %>% n_distinct() #only 259/512 (-259 not seeded or found) have BTRE
# allwy23 <- allwy23 %>% mutate(log.inv = log(invcov)) %>% 
#   mutate(log.inv = ifelse(log.inv == -Inf, NA, log.inv))

#### remove all plots where CWM could not be validly calculated
suballwy23 <- allwy23 %>% filter(propnativebrte >= 80)

#model
summary(lmer(log.propbrte~trt*drought+ (1|block), allwy23))
summary(irmod <- lmer(log.propbrte~distir*drought+ (1|block), suballwy23))
summary(m2<-lmer(log.propbrte~distdt*drought+ (1|block), suballwy23))
summary(lmer(log.propbrte~distfd*drought+ (1|block), allwy23))

# summary(irmod <- glmer(BRTEpres~trt*drought+ (1|block), fortest, family = "binomial"))
# x <- ggpredict(irmod,c("trt","drought")) 
# plot(x, show.title = F)

anova(lmer(log.gr~distdt*trt*distfd*drought*year+ (1|block), testno))

invboxwy <- ggplot(allwy23, aes(y=log.propbrte,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  #facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y=" ", fill="seed trt")+
  theme_ggeffects()
suballwy23 <- suballwy23 %>% mutate(nativecovbin = cut(nativecov.plot, 
                                                       breaks = quantile(nativecov.plot, probs = seq(0, 1, length.out = 4), na.rm = TRUE), 
                                                       include.lowest = TRUE))
# mutate(nativecovbin = bins.quantiles(nativecov.plot, 3, 3))
#   nativecovbin = ifelse(nativecov.plot 
invirwy <- ggplot(suballwy23, aes(y=log.propbrte,x=distir,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  #facet_wrap(~year)+
  labs(y=" ", x="IR target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
invdtwy <- ggplot(suballwy23, aes(y=log.propbrte,x=distdt,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  #facet_wrap(~year)+
  labs(y=" ", x="DT target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
invfdwy <- ggplot(suballwy23, aes(y=log.propbrte,x=distfd,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  #facet_wrap(~year)+
  labs(y=" ", x="FD target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distrwy <- ggplot(suballwy23, aes(y=log.brte,x=distr,col=drought))+
  geom_point(aes(y=log.brte,x=distr,col=drought, shape=year))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  #facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(geom = "label",label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_te
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()


#### combined figures
library(ggpubr)
wyfigtop <- ggarrange(invboxwy,invirwy, 
                      common.legend = T, legend = "right",
                      labels = c("a","b"),label.x = 1)
wyfigbottom <-ggarrange(invfdwy,invdtwy, 
                        common.legend = T, legend = "right",
                        labels = c("c","d"),label.x = 1)
wyfiginvasion <- ggarrange(wyfigtop,wyfigbottom, nrow=2)
wyfiginvasion <- annotate_figure(wyfiginvasion, bottom = "Euclidean distance to CWM target",
                                 left="log(relative cover BRTE)")

tiff("figures/invasionfigwy.tiff", res=400, height = 5,width =8, "in",compression = "lzw")
wyfiginvasion
dev.off()