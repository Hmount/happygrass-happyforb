#### Log-response ratios

library(tidyverse)


# start with CA
## load in data, clean and modify columns
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
droughtcols <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

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

hist(test$log.gr)
# 
# testno20 <- test %>% filter(!is.na(log.gr))
testno <- test %>% filter(year!="2021")
testonly22 <- test %>% filter(year=="2022")
m3.ca <- lmer(log.gr ~ trt * distdt * drought * year + (1|structure), data = testno)
summary(m3.ca)
anova(m3.ca) #cov effected by drought, and dist:drought int.
emm.ca <- emmeans(m3.ca, c("distdt","trt","drought"), at = list(distdt = 0.56))#,"dist"))
pairs(emm.ca)#,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.ca,pairwise~"drought","trt")
contrast(emm.ca)

difflsmeans(m3.ca, at = list(distdt = 0.56))


#make pairwise comparisons and letters to display
emm <- emmeans(m3.ca, specs = ~ trt * drought * distdt * year,  at = list(distdt = c(0.56,1.56)))
pairwise <- pairs(emm)
summary(pairwise, adjust = "tukey")
letters <- multcomp::cld(emm, Letters = LETTERS)

summary(t <- aov(distdt~trt*drought*year, suballwy))
tuktest <- TukeyHSD(t)

letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought:year'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\2", rownames(letterstest)))
letterstest$year <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\3", rownames(letterstest)))
test <- allwy %>% group_by(year, drought, trt) %>% summarise(yposition = quantile(dist,.8))
test <- merge(letterstest,test, by = c("year", "drought", "trt"))
test2 <- merge(test,allwy, by = c("year", "drought", "trt"), all=T)


#multcompView::multcompLetters4(t,tuktest)
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- allwy %>% group_by(drought, trt) %>% summarise(yposition = quantile(distdt,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,allwy, by = c("drought", "trt"), all=T)

ggplot(test, aes(y=log.gr,x=trt,fill=drought))+
  #geom_point()+
  #geom_point(alpha=.3, pch=20)+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)+
  facet_wrap(~year,scales="fixed")+
  # geom_text(aes(y=yposition,label = Letters), 
  #           position = position_dodge(width = 0.9), 
  #           vjust = -0.5,
  #           #angle = 15,
  #           size=3) +
  theme_classic()

gr.raw.fig.ca <- ggplot(testno, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y="(log) annual growth rate", fill="seed trt")+
  theme_ggeffects()
x <- ggpredict(m3.ca,c("trt","distdt [c(.56,1.56)]","drought","year")) 
gr.pred.fig.ca <- plot(x, show.title = F)+ labs(col="growth rate")

# x <- ggpredict(m3.ca,c("year","distdt [c(.56,1.56)]","trt","drought")) 
# plot(x, show.title = F, connect_lines = T, one_plot = T)


#rootdiam model
#suballca$rootdiam <- normalize(suballca$rootdiam) # normalize FD
m2rd.ca2 <- lmer(log.gr ~ drought * distdt * Year+ (1 | structure), data = test)
m2rd.ca2 <- lmer(log.gr ~ distdt + (1|year)+(1 | structure), data = test)
summary(m2rd.ca2) 
anova(m2rd.ca2) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
testplot <- test %>% filter(!is.na(log.gr))
gr.dist.fig.ca <- ggplot(testplot, aes(y=log.gr,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
ggplot(test, aes(y=log.gr,x=distdt))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(~Year)

#together
ggarrange(gr.raw.fig.ca,gr.dist.fig.ca,gr.pred.fig.ca, nrow=3)

#figure 2/8
distmod <- lmer(log.gr ~ distir * year + (1 | block), data = testno)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(testno, aes(y=log.gr,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  #stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  labs(y="absolute cover native")+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
dtca<-ggplot(testno, aes(y=log.gr,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  stat_cor(label.x = 1.6)+
  # labs(y="relative cover BRTE")+
  theme_ggeffects()
irca<-ggplot(testno, aes(y=log.gr,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  #labs(y="relative cover BRTE")+
  facet_wrap(~year)+
  stat_cor()+
  theme_ggeffects()+
  theme(legend.position = "none")
fdca<-ggplot(testno, aes(y=log.gr,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  #labs(y="relative cover BRTE")+
  facet_wrap(~year)+  stat_cor(label.y=c(4,3))+
  theme_ggeffects()+
  theme(legend.position = "none")
#distplotswy<-ggarrange(irwy,dtwy,fdwy, nrow=2, ncol=2,common.legend = T )

#plots 2/8
boxca <-ggplot(testno, aes(y=log.gr ,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="annual growth rate", fill="seed trt")+
  theme_ggeffects()

ggarrange(irca,dtca,fdca, boxca,nrow=2, ncol=2,common.legend = F)











#### now WY
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
#wydist <- wydist %>% select(-X) #%>% filter(trt.b.sub.y!="target")
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

hist(test$log.gr)

testno20 <- test %>% filter(!is.na(log.gr))
testno <- test %>% filter(year!="2021")
testonly22 <- test %>% filter(year=="2023")
m3.wy <- lmer(log.gr ~ trt * distdt * drought + (1|block), data = testonly22)
summary(m3.wy)
anova(m3.wy)
emm.ca <- emmeans(m3.ca, c("distdt","trt","drought"), at = list(distdt = 0.56))#,"dist"))
pairs(emm.ca)#,simple="drought") #at an average distance from target how does cover vary as a result of the other vars
emmeans(m3.ca,pairwise~"drought","trt")
contrast(emm.ca)

difflsmeans(m3.ca, at = list(distdt = 0.58))


#make pairwise comparisons and letters to display
emm <- emmeans(m3.wy, specs = ~ trt * drought * distdt,  at = list(distdt = c(0.56,1.56)))
pairwise <- pairs(emm)
summary(pairwise, adjust = "tukey")
letters <- multcomp::cld(emm, Letters = LETTERS)

ggplot(testonly22, aes(y=log.gr,x=trt,fill=drought))+
  #geom_point()+
  #geom_point(alpha=.3, pch=20)+
  geom_boxplot()+
  scale_fill_manual(values=droughtcols)+
 # facet_wrap(~year,scales="free")+
  theme_classic()

x <- ggpredict(m3.wy,c("trt","distdt [c(.48,2.84)]","drought")) 
plot(x)#+
  labs(x="seeding treatment", 
       y="absolute cover native species",
       title=" ")#+
  #scale_color_manual(values = c("1" = "slateblue", "2" = "lightsalmon3"), #orchid4 #olivdrab3
                     labels = c("1" = "0.48", "2" = "2.84"))#+
#geom_text(data = letters, aes(x = trt, color=distdt, y = emmean, label = .group), hjust=2, vjust=-1.75, col="black") #add tukey labels


#### 2/7/ comparison to loggr
testno <- test %>% filter(year!="2021")
#testonly22 <- test %>% filter(year=="2022")
m3.wy <- lmer(log.gr ~ trt * distdt * drought * year + (1|block), data = testno)
summary(m3.wy)
anova(m3.wy) 

dissboxwy <- ggplot(testno, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y="log(annual growth rate)", fill="seed trt")+
  theme_ggeffects()
x <- ggpredict(m3.wy,c("trt","distdt [c(.48,2.84)]","drought", "year")) 
plot(x, show.title = F)

distmod <- lmer(log.gr ~ distir * year + (1 | block), data = testno)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
wyplotdiss <- ggplot(testno, aes(y=log.gr,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  #stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  labs(y="log(annual growth rate)")+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

#figure 2/8 WY
distmod <- lmer(log.gr ~ distir * year + (1 | block), data = testno)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
ggplot(testno, aes(y=log.gr,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  #stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  labs(y="absolute cover native")+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
dtwy<-ggplot(testno, aes(y=log.gr,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  stat_cor()+#label.x = 1.85)+
 # labs(y="relative cover BRTE")+
  theme_ggeffects()
irwy<-ggplot(testno, aes(y=log.gr,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  #labs(y="relative cover BRTE")+
  facet_wrap(~year)+
  stat_cor()+
  theme_ggeffects()+
  theme(legend.position = "none")
fdwy<-ggplot(testno, aes(y=log.gr,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  #labs(y="relative cover BRTE")+
  facet_wrap(~year)+  stat_cor(label.y=c(4,3))+
  theme_ggeffects()+
  theme(legend.position = "none")
#distplotswy<-ggarrange(irwy,dtwy,fdwy, nrow=2, ncol=2,common.legend = T )

#plots 2/8
boxwy <-ggplot(testno, aes(y=log.gr ,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="annual growth rate", fill="seed trt")+
  theme_ggeffects()

ggarrange(irwy,dtwy,fdwy, boxwy,nrow=2, ncol=2,common.legend = F)


######### exploring escape stategy
modcwmca <- lmer(log.gr ~ drought * graminoid * Year+ (1 | structure), data = testno)
anova(modcwmca)
ggplot(testno, aes(y=log.gr,x=graminoid,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  theme_ggeffects()
#n and lma alright model, maybe rd, SRL almost all, RMF and seedmass no 



modcwmwy <- lmer(log.gr ~ drought * graminoid * year+ (1 | block), data = testno)
anova(modcwmwy)
ggplot(testno, aes(y=log.gr,x=graminoid,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  facet_wrap(~year)+
  theme_ggeffects()
#n and lma alright model, maybe rd, SRL almost all, RMF and seedmass no 
















#### #all combined 3/7/24 for dissertation proposal ####
# start with CA
## load in data, clean and modify columns
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
droughtcols <- c("cntl"="skyblue", "drt"="tomato") #create variable for color

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
test22 <- test %>% filter(year=="2022")

#model
anova(lmer(log.gr~distdt*trt*drought+(1|year)+ (1|block), testno))
anova(lmer(log.gr~distdt*trt*drought*year+ (1|block), testno))

dissboxca <- ggplot(testno, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  #facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y="", fill="seed trt")+
  theme_ggeffects()
testplot <- testno %>% filter(!is.na(log.gr))
caplotdiss <- ggplot(testplot, aes(y=log.gr,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y=" ", x=" ")+
  facet_wrap(~year)+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  geom_hline(yintercept =0,col="black")+
  theme_ggeffects()



#### now WY
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
testnoD <- testno %>% filter(drought=="cntl")

#model
summary(lmer(log.gr~distdt*trt*drought*year+ (1|block), testno))
anova(lmer(log.gr~distdt*trt*year+ (1|block), testnoD))

dissboxwy <- ggplot(testnoD, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  facet_wrap(~year,scales="fixed")+
  geom_hline(yintercept =0,col="black")+
  labs(y=" ", fill="seed trt")+
  theme_ggeffects()
wyplotdiss <- ggplot(testnoD, aes(y=log.gr,x=distdt,color=year))+
  geom_point()+
  geom_smooth(method = "lm")+
  #scale_color_manual(values=droughtcols)+
  #facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()


library(ggpubr)
dissfig1 <-ggarrange(dissboxwy,dissboxca, common.legend = T, legend = "right")
dissfig2 <-ggarrange(wyplotdiss, caplotdiss, common.legend = T,legend = "right")
dissfig2 <- annotate_figure(dissfig2, bottom = "Euclidean distance to drought tolerant CWM targets")
dissfig <-ggarrange(dissfig1,dissfig2, nrow=2)
testfig <- annotate_figure(dissfig, left = "log(annual growth rate)", top = "WY                                                    CA")

