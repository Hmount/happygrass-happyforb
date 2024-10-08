#### Analysis of invasion, CA
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
#cadist <- cadist %>% select(-X) #%>% filter(trt.b.y!="target")
#break apart distances ID to make wider and merge together
cadist <- separate(cadist, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist <- cadist %>% filter(trt!="target")

#cadist <- cadist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allca <- merge(allca,cadist, 
               by=c("year","trt","block"), all=T)
allca23 <- allca %>% filter(year=="2023") #2023 makes sense to use because all grasses had equal seeding
allca23$trt <- as.factor(allca23$trt)
allca23$plot.y <- as.factor(allca23$plot.y)
allca23year <- as.factor(allca23$year)

### how many WY 2022+2023 data need to be dropped from CWM calculations 
allca23 %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  filter(propnative < 80)
18/204*100 # only 9% total observation to remove for inv. models
allca23 <- allca23 %>% 
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  mutate(propfesper = FESPER/(native.cover+inv.grass.cov)) %>% #this one is not converted to percent
  mutate(propnativefesper = propnative+(FESPER*100))

#make dought column
allca23 <- allca23 %>% mutate(drought = as.factor(ifelse(water=="0.5","drt","cntl")))

## Ensure levels are correctly compared in models
allca23$trt <- relevel(allca23$trt, ref = "rand") #make random communities the reference level
allca23$drought <- relevel(allca23$drought, ref = "cntl") #make ambient precip the reference level
#for visuals
droughtcolsca <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

### data summary
# look at response variable in each dataset
# CA
hist(allca23$FESPER)
hist(log(allca23$FESPER)) #better logged
allca23 %>% filter(FESPER!="0") %>% n_distinct() #72/204 (-100 not seeded or found) have FESPER
allca23 <- allca23 %>% mutate(log.fesper = log(FESPER)) %>%  mutate(log.fesper = ifelse(log.fesper == -Inf, NA, log.fesper))
hist(allca23$propfesper)
hist(log(allca23$propfesper)) #better logged
allca23 %>% filter(propfesper!="0") %>% n_distinct() #72/204 (-100 not seeded or found) have FESPER
allca23 <- allca23 %>% mutate(log.propfesper = log(propfesper)) %>%  mutate(log.propfesper = ifelse(log.propfesper == -Inf, NA, log.propfesper))
# # I am not convinced of this variable, but has more even spread across precipitation treatments 
hist(allca23$inv.grass.cov)
hist(log(allca23$inv.grass.cov)) #better logged
allca23 %>% filter(inv.grass.cov!="0") %>% n_distinct() #only 87/204 (-100 not seeded or found) have FESPER
allca23 <- allca23 %>% mutate(log.invg = log(inv.grass.cov)) %>%
  mutate(log.invg = ifelse(log.invg == -Inf, NA, log.invg))

# assess only invaded plots
suballca23 <- allca23 %>% filter(fesper.seeded=="1")
suballca23 <- suballca23 %>% filter(mono=="0") #no mono

#### remove all plots where CWM could not be validly calculated
suballca23x <- suballca23 %>% filter(propnativefesper >= 80)


#model
summary(trtmod<-lmer(log.invg~trt*drought+ (1|block), suballca23))
summary(irmod <- lmer(log.invg~distir*drought+ (1|block), suballca23))
summary(fdmod<-lmer(log.invg~distfd*drought+ (1|block), suballca23))
summary(dtmod<-lmer(log.invg~distdt*drought+ (1|block), suballca23))

# summary(irmod <- glmer(BRTEpres~trt*drought+ (1|block), fortest, family = "binomial"))
# x <- ggpredict(irmod,c("trt","drought")) 
# plot(x, show.title = F)

anova(lmer(log.gr~distdt*trt*distfd*drought*year+ (1|block), testno))

invboxca <- ggplot(suballca23, aes(y=log.invg ,x=drought,fill=trt))+
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
invirca <- ggplot(suballca23, aes(y=log.invg,x=distir,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca)+
 # facet_wrap(~year)+
  labs(y=" ", x="IR target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
invdtca <- ggplot(suballca23, aes(y=log.invg,x=distdt,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  #facet_wrap(~year)+
  labs(y=" ", x="DT target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
invfdca <- ggplot(suballca23, aes(y=log.invg,x=distfd,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  #facet_wrap(~year)+
  labs(y=" ", x="FD target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()


#### combined figures
library(ggpubr)
cafigtop <- ggarrange(invboxca,invirca, 
                      common.legend = T, legend = "right",
                      labels = c("a","b"),label.x = 1)
cafigbottom <-ggarrange(invfdca,invdtca, 
                        common.legend = T, legend = "right",
                        labels = c("c","d"),label.x = 1)
cafiginvasion <- ggarrange(cafigtop,cafigbottom, nrow=2)
cafiginvasion <- annotate_figure(cafiginvasion, bottom = "Euclidean distance to CWM target",
                                left="log(relative cover invasive grass)")

tiff("figures/invasionfigca.tiff", res=400, height = 5,width =8, "in",compression = "lzw")
cafiginvasion
dev.off()