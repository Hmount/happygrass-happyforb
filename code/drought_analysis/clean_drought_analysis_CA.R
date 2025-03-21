#### Analysis of drought treatments, CA
#### how does drought effect native species in the different seeded communities? and
#### with our distance from targets?
#### Are DT communities more tolerant of drought than random? than FD? 
#### Are communities with more DT-target traits more tolerant of drought?

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
comp.ca$structure <- as.factor(comp.ca$structure)
comp.ca <- comp.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
#save comp.ca$struture to add onto other datframes and use as spatial variable
ca.structure <- comp.ca %>% select(c("block","trt","structure"))
ca.structure <- distinct(ca.structure)

cwm.ca <- read.csv("data/cwm_ca.csv")# California CWM data
cwm.ca <- cwm.ca %>% filter(year != "0") #remove predicted/target communities
cwm.ca$year <- as.factor(cwm.ca$year)
cwm.ca$block <- as.factor(cwm.ca$block)
#make new sequence column
cwm.ca <- cwm.ca %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0"))))
cwm.ca$yrorder <- as.numeric(cwm.ca$yrorder)
cwm.ca <- merge(cwm.ca, ca.structure, by=c("block", "trt"), all.x = T)
#add plot ID column (but give NA to target/predicted communities)
cwm.ca <- cwm.ca %>% 
  mutate(plot = paste(block, trt, sep = "."))
cwm.ca <- cwm.ca %>% filter(year != "0") #remove predicted/target communities

cwm.ca <- cwm.ca %>% select(-Rdiam) #remove CWM rootdiam column to avoid confusion
cwmFD <- read.csv("data/cwm_raoq_ca.csv") #add FD for traits that need it (rootdiam)
cwmFD <- cwmFD %>% select(block,trt,year,water,rootdiam, full) #only columns we need
cwmFD <- merge(cwmFD, ca.structure, all.x = T)
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

#make precipitation treatment column
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
cadat <- merge(suballca, forgrate, all.x=T)

#find annual growth rate
cadat <- cadat %>% mutate(growrate = native.cover/covprevyr)
cadat$log.gr <- log(cadat$growrate) 
cadatno21 <- cadat %>% filter(year!="2021")
#test22 <- test %>% filter(year=="2022")

#model
summary(trtmod <- lmer(log.gr~trt*drought*year+ (1|structure.x), cadatno21))
summary(dtmod<-lmer(log.gr~distdt*drought*year+ (1|structure.x), cadatno21))
summary(fdmod<-lmer(log.gr~distfd*drought*year+ (1|structure.x), cadatno21))
summary(irmod<-lmer(log.gr~distir*drought*year+ (1|structure.x), cadatno21))

#### Making plots
## create label names for facets used in all plots:
labelnames.ca <- c('2022' = "2021 - 2022",
                   '2023' = "2022 - 2023")

## ~ seeding treatment
#create letters for plotting:
library(emmeans)
# Step 1: Get the emmeans for the interaction of trt, drought, and year
emm_trt <- emmeans(trtmod, ~ trt * drought * year)
## # Step 2: Obtain pairwise contrasts for the interaction
## contrast_trt <- contrast(emm_trt, method = "pairwise")
# Step 2: Generate the compact letter display using multcomp::cld
letters <- multcomp::cld(emm_trt, alpha = 0.05, Letters = letters, adjust = "tukey")
# Step 3: Convert the results to a data frame
letters_df <- as.data.frame(letters)
# Step 4: Create a temporary data frame with the desired y-position for plotting
dttemp2 <- cadatno21 %>%
  group_by(drought, trt, year) %>%
  summarise(yposition = max(log.gr, na.rm=T), .groups = 'drop')
# Step 5: Merge the letter results with the y-position data
dttemp2 <- merge(letters_df, dttemp2, by = c("drought", "trt","year"))
# Merge with the original data to get the final dataset
dttemp3 <- merge(cadatno21, dttemp2, by = c("drought", "trt", "year"), all = TRUE)

#plot:
dissboxca <- ggplot(dttemp3, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  geom_text(aes(y=yposition,label = .group), 
            position = position_dodge(width = .9), 
            vjust = -0.5,
            size=3)+
  ylim(c(NA,3.5))+
  #scale_fill_manual(values = c("#482576B3","#2A788EB3","#43BF71B3"))+ #viridis::viridis(4, option="D",begin = .1, end = 1, alpha = 0.7)
  scale_x_discrete(labels = c("Addition", "Reduction"))+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7,
                      labels = c("RC","DT","FD","IR"))+
  facet_wrap(~year,scales="fixed", labeller = as_labeller(labelnames.ca))+
  geom_hline(yintercept =0,col="black")+
  labs(y=" ", fill="Seeding 
Treatment", x = "Precipitation treatment")+
  theme_ggeffects()

## ~ distance to DT traits
distdtca <- ggplot(cadatno21, aes(y=log.gr,x=distdt,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca, labels = c("Addition", "Reduction"))+
  facet_wrap(~year, labeller = as_labeller(labelnames.ca))+
  labs(y=" ", x="Euclidean distance to DT target", col="Precipitation 
Treatment")+
  ylim(c(NA,3.5))+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

## ~ distance to IR traits
distirca <- ggplot(cadatno21, aes(y=log.gr,x=distir,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca, labels = c("Addition", "Reduction"))+
  facet_wrap(~year, labeller = as_labeller(labelnames.ca))+
  labs(y=" ", x="Euclidean distance to IR target", col="Precipitation 
Treatment")+
  ylim(c(NA,3.5))+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

## ~ distance to FD traits
distfdca <- ggplot(cadatno21, aes(y=log.gr,x=distfd,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolsca, labels = c("Addition", "Reduction"))+
  facet_wrap(~year, labeller = as_labeller(labelnames.ca))+
  labs(y=" ", x="Euclidean distance to FD target", col="Precipitation 
Treatment")+
  ylim(c(NA,3.5))+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(2.5,2.4),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

#### combined figures
library(ggpubr)
cafigleft <- ggarrange(dissboxca,distfdca, nrow=2,
                      common.legend = T, legend = "bottom",
                      labels = c("a","c"),label.x = .05)
cafigright <-ggarrange(distdtca,distirca, nrow=2,
                        common.legend = T, legend = "bottom",
                        labels = c("b","d"),label.x = .05)
cafigdrought <- ggarrange(cafigleft,cafigright, ncol=2)
cafigdrought <- annotate_figure(cafigdrought, 
                                left="Annual growth rate")

tiff("figures/droughtfigca.tiff", res=400, height = 5,width =8, "in",compression = "lzw")
cafigdrought
dev.off()