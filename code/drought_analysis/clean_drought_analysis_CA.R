#### Analysis of drought treatments, CA
#### how does drought effect native species in the different seeded communities? and
#### with our distance from targets?
#### Are DT communities more tolerant of drought than random? than FD? 
#### Are communities with more DT-target traits more tolerant of drought?

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## read in all data and clean to create master dataframe
comp.ca <- read.csv("data/comp_ca.csv") #read in California comp data
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
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

cwm.ca <- cwm.ca %>% select(-Rdiam) #remove CWM rootdiam column to avoid confusion
cwmFD <- read.csv("data/cwm_raoq_ca.csv") #add FD for traits that need it (rootdiam)
cwmFD <- cwmFD %>% select(block,trt,year,water,rootdiam, full) #only columns we need
cwm.ca <- merge(cwm.ca,cwmFD, all.x=T)

# combine to master df (remove spp columns for now)
allca <- merge(comp.ca[,-c(18:53,55:66)],cwm.ca,  #this merge drops monoculture plots
               by=c("year","trt","block","water"))#, 
allca$trt <- as.factor(allca$trt)

# also combine CWM_distances dataframe to master df 
cadist <- read.csv("data/cwm_maxdistances_ca.csv")
#break apart distances ID to make wider and merge together
cadist <- separate(cadist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist <- cadist %>% filter(trt!="target")

allca <- merge(allca,cadist, 
               by.y=c("year","trt","block"), all=T)
allca$trt <- as.factor(allca$trt)
allca$plot.y <- as.factor(allca$plot.y)
allca$Year <- as.factor(allca$year)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
test <- allca %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  filter(propnative < 80)
table(test$year)
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
droughtcolsca <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# allca <- allca %>%
#   mutate(nativecov_tran = sqrt(native.cover))

#### remove all plots where CWM could not be validly calculated
suballca <- allca %>% filter(propnative >= 80)


#attach previous years cover as a column to calculate response ratio  
grate21 <- suballca %>% filter(year=="2021")
grate21 <- grate21 %>% select(c(block,trt,native.cover))
colnames(grate21) <- c("block", "trt","covprevyr")
grate21$year <- "2022"
grate22 <- suballca %>% filter(year=="2022")
grate22 <- grate22 %>% select(c(block,trt,native.cover))
colnames(grate22) <- c("block", "trt","covprevyr")
grate22$year <- "2023"
forgrate <- bind_rows(grate21,grate22)
test <- merge(suballca, forgrate, all.x=T)

#find annual growth rate
test <- test %>% mutate(growrate = native.cover/covprevyr)
test$log.gr <- log(test$growrate) 
testno <- test %>% filter(year!="2021")
#test22 <- test %>% filter(year=="2022")

#model
anova(lmer(log.gr~distdt*trt*drought+(1|year)+ (1|block), testno))
anova(lmer(log.gr~distdt*trt*drought*year+ (1|block), testno))

dissboxca <- ggplot(testno, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y="", fill="seed trt")+
  theme_ggeffects()
testplot <- testno %>% filter(!is.na(log.gr))
dissdtca <- ggplot(testno, aes(y=log.gr,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca)+
  labs(y=" ", x=" ")+
  facet_wrap(~year)+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
dissirca <- ggplot(testno, aes(y=log.gr,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca)+
  labs(y=" ", x=" ")+
  facet_wrap(~year)+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
dissfdca <- ggplot(testno, aes(y=log.gr,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca)+
  labs(y=" ", x=" ")+
  facet_wrap(~year)+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
dissrca <- ggplot(testno, aes(y=log.gr,x=distr,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca)+
  labs(y=" ", x=" ")+
  facet_wrap(~year)+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
