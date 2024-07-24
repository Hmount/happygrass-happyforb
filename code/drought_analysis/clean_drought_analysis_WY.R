#### Analysis of drought treatments, WY 
#### how does drought effect native species in the different seeded communities? and
#### with our distance from targets?
#### Are DT communities more tolerant of drought than random? than FD? 
#### Are communities with more DT-target traits more tolerant of drought?

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## read in all data and clean to create master dataframe
comp.wy <- read.csv("data/comp_wy_plot.csv") # Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
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
allwy <- merge(comp.wy[,-c(5:60)],cwm.wy, by=c("year","trt","block"))
allwy$trt <- as.factor(allwy$trt)

# also combine CWM_distances dataframe to master df 
wydist <- read.csv("data/cwm_maxdistances_wy(plot).csv")
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
#for visuals
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

#### remove all plots where CWM could not be validly calculated
suballwy <- allwy %>% filter(propnative >= 80)

#attach previous years cover as a column to calculate response ratio  
grate21 <- suballwy %>% filter(year=="2021")
grate21 <- grate21 %>% select(c(block,trt,nativecov.plot))
colnames(grate21) <- c("block", "trt","covprevyr")
grate21$year <- "2022"
grate22 <- suballwy %>% filter(year=="2022")
grate22 <- grate22 %>% select(c(block,trt,nativecov.plot))
colnames(grate22) <- c("block", "trt","covprevyr")
grate22$year <- "2023"
forgrate <- bind_rows(grate21,grate22)
test <- merge(suballwy, forgrate, all.x=T)

#find annual growth rate
test <- test %>% mutate(growrate = nativecov.plot/covprevyr)
test$log.gr <- log(test$growrate) 
testno <- test %>% filter(year!="2021")
#testnoD <- testno %>% filter(drought=="cntl")

#model
summary(lmer(log.gr~distdt*trt*drought*year+ (1|block), testno))
anova(lmer(log.gr~distdt*trt*drought*year+ (1|block), testno))

dissboxwy <- ggplot(testno, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y=" ", fill="seed trt")+
  theme_ggeffects()
distdtwy <- ggplot(testno, aes(y=log.gr,x=distdt,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distirwy <- ggplot(testno, aes(y=log.gr,x=distir,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distfdwy <- ggplot(testno, aes(y=log.gr,x=distfd,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distrwy <- ggplot(testno, aes(y=log.gr,x=distr,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()


hist(test$log.gr)

difflsmeans(m3.ca, at = list(distdt = 0.58))

library(ggpubr)
dissfig1 <-ggarrange(dissboxwy,dissboxca, common.legend = T, legend = "right")
dissfig2 <-ggarrange(wyplotdiss, caplotdiss, common.legend = T,legend = "right")
dissfig2 <- annotate_figure(dissfig2, bottom = "Euclidean distance to drought tolerant CWM targets")
dissfig <-ggarrange(dissfig1,dissfig2, nrow=2)
testfig <- annotate_figure(dissfig, left = "log(annual growth rate)", top = "WY                                                    CA")
################################################################

