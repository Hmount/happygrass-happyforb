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
anova(lmer(log.gr~trt*drought*year+ (1|block), testno))
summary(dtmod <- lmer(log.gr~distdt*drought*year+ (1|block), testno))
summary(lmer(log.gr~distfd*drought*year+ (1|block), testno))
summary(lmer(log.gr~distir*drought*year+ (1|block), testno))

anova(lmer(log.gr~distdt*trt*distfd*drought*year+ (1|block), testno))

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
  scale_color_manual(values=droughtcolswy)+
  facet_wrap(~year)+
  labs(y=" ", x="DT target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distirwy <- ggplot(testno, aes(y=log.gr,x=distir,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  facet_wrap(~year)+
  labs(y=" ", x="IR target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distfdwy <- ggplot(testno, aes(y=log.gr,x=distfd,col=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  facet_wrap(~year)+
  labs(y=" ", x="FD target")+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
distrwy <- ggplot(testno, aes(y=log.gr,x=distr,col=drought))+
  geom_point(aes(y=log.gr,x=distr,col=drought, shape=year))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy)+
  facet_wrap(~year)+
  labs(y=" ", x="")+
  geom_hline(yintercept =0,col="black")+
  stat_cor(geom = "label",label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_te
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()

hist(test$log.gr)

difflsmeans(m3.ca, at = list(distdt = 0.58))


#### combined figures
wyfigtop <- ggarrange(dissboxwy,distdtwy, 
                      common.legend = T, legend = "right",
                      labels = c("a","b"),label.x = 1)
wyfigbottom <-ggarrange(distfdwy,distirwy, 
                        common.legend = T, legend = "right",
                        labels = c("c","d"),label.x = 1)
wyfigdrought <- ggarrange(wyfigtop,wyfigbottom, nrow=2)
wyfigdrought <- annotate_figure(wyfigdrought, bottom = "Euclidean distance to CWM target",
                                left="Annual growth rate (log(lambda))")


tiff("figures/droughtfigwy.tiff", res=400, height = 5,width =8, "in",compression = "lzw")
wyfigdrought
dev.off()

tiff("figures/ESAfig3.tiff", res=400, height = 6,width =8, "in",compression = "lzw")
wyfigdrought
dev.off()

dissfig1 <-ggarrange(dissboxwy,distdtwy,distfdwy,distirwy, 
                     common.legend = F, legend = "right")
dissfig2 <- annotate_figure(dissfig1, bottom = "Euclidean distance to drought tolerant CWM targets")

#### for ESA pres
ESAdissfig1 <-ggarrange(dissboxwy,distdtwy, 
                     common.legend = T, legend = "right")
ESAdissfig2 <- annotate_figure(ESAdissfig1, bottom = "Euclidean distance to drought tolerant CWM targets")




library(ggpubr)
dissfig1 <-ggarrange(dissboxwy,dissboxca, common.legend = T, legend = "right")
dissfig2 <-ggarrange(wyplotdiss, caplotdiss, common.legend = T,legend = "right")
dissfig2 <- annotate_figure(dissfig2, bottom = "Euclidean distance to drought tolerant CWM targets")
dissfig <-ggarrange(dissfig1,dissfig2, nrow=2)
testfig <- annotate_figure(dissfig, left = "log(annual growth rate)", top = "WY                                                    CA")
################################################################

#####output of models as formatted tables:
library(lme4)
library(broom.mixed)
library(knitr)
library(kableExtra)

model_summary <- tidy(dtmod)

# Separate fixed and random effects
fixed_effects <- model_summary %>% filter(effect == "fixed")
random_effects <- model_summary %>% filter(effect == "ran_pars")

# Adjusting Terms to be More Descriptive
fixed_effects <- fixed_effects %>%
  mutate(term = recode(term,
                       `(Intercept)` = "Intercept",
                       `trt1` = "Treatment 1",
                       `trt2` = "Treatment 2",
                       `trt3` = "Treatment 3",
                       `drought1` = "Drought Level 1",
                       `drought2` = "Drought Level 2",
                       `year1` = "Year 1",
                       `year2` = "Year 2"))

random_effects <- random_effects %>%
  mutate(term = recode(term,
                       `(Intercept)` = "Intercept"))

# Combine fixed and random effects with proper column names
fixed_effects <- fixed_effects %>%
  select(term, estimate, std.error, statistic, p.value) %>%
  mutate(effect_type = "Fixed Effect")

random_effects <- random_effects %>%
  select(term, estimate, std.error) %>%
  mutate(statistic = NA, p.value = NA, effect_type = "Random Effect")

# Rename columns to be more descriptive
colnames(fixed_effects) <- c("Term", "Estimate", "Std. Error", "Statistic", "P-value", "Effect Type")
colnames(random_effects) <- c("Term", "Estimate", "Std. Error", "Statistic", "P-value", "Effect Type")

# Combine both into one dataframe
combined_effects <- bind_rows(fixed_effects, random_effects)

# Create the table
table<-kable(combined_effects, caption = "Summary of Linear Mixed-Effects Model", digits = 3) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"),full_width = F, position = "left")

# Save the table as HTML
save_kable(table, "model_summaryirca.doc")
