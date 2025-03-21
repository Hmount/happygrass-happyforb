#### Analysis of drought treatments, WY 
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
wydat <- merge(suballwy, forgrate, all.x=T)

#find annual growth rate
wydat <- wydat %>% mutate(growrate = nativecov.plot/covprevyr)
wydat$log.gr <- log(wydat$growrate) 
wydatno21 <- wydat %>% filter(year!="2021")
#testnoD <- testno %>% filter(drought=="cntl")

#model
summary(trtmod<-lmer(log.gr~trt*drought*year+ (1|block), wydatno21))
summary(dtmod <- lmer(log.gr~distdt*drought*year+ (1|block), wydatno21))
summary(fdmod <- lmer(log.gr~distfd*drought*year+ (1|block), wydatno21))
summary(irmod <- lmer(log.gr~distir*drought*year+ (1|block), wydatno21))

#anova(lmer(log.gr~distdt*trt*distfd*drought*year+ (1|block), testno))


#### Making plots
## create label names for facets used in all plots:
labelnames.wy <- c('2022' = "2021 - 2022 (dry)",
                   '2023' = "2022 - 2023 (wet)")

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
dttemp2 <- wydatno21 %>%
  group_by(drought, trt, year) %>%
  summarise(yposition = quantile(log.gr,.8, na.rm = T), .groups = 'drop')
# Step 5: Merge the letter results with the y-position data
dttemp2 <- merge(letters_df, dttemp2, by = c("drought", "trt","year"))
# Merge with the original data to get the final dataset
dttemp3 <- merge(wydatno21, dttemp2, by = c("drought", "trt", "year"), all = TRUE)

#plot:
dissboxwy <- ggplot(dttemp3, aes(y=log.gr,x=drought,fill=trt))+
  geom_boxplot()+
  geom_text(aes(y=yposition,label = .group), 
            position = position_dodge(width = .85), 
            vjust = -1.75,
            hjust = .6,
            size=3)+
  scale_x_discrete(labels = c("Ambient", "Reduction"))+
  ylim(c(-5,5))+
  #scale_fill_manual(values = c("#482576B3","#2A788EB3","#43BF71B3"),
  #                  labels = c("RC","DT","FD"))+ #viridis::viridis(4, option="D",begin = .1, end = 1, alpha = 0.7)
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7,
                      labels = c("RC","DT","FD","IR"))+
  facet_wrap(~year,scales="fixed", labeller = as_labeller(labelnames.wy))+
  geom_hline(yintercept =0,col="black")+
  labs(y=" ", fill="Seeding 
Treatment", x = "Precipitation treatment")+
  theme_ggeffects()


## ~ distance to DT traits
distdtwy <- ggplot(wydatno21, aes(y=log.gr,x=distdt,col=drought))+
  geom_point(aes(alpha=.8))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy, labels = c("Ambient", "Reduction"))+
  facet_wrap(~year, labeller = as_labeller(labelnames.wy))+
  labs(y=" ", x="Euclidean distance to DT target", col="Precipitation 
Treatment")+
  ylim(c(-5,5))+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  guides(alpha = "none")+
  theme_ggeffects()

## ~ distance to IR traits
distirwy <- ggplot(wydatno21, aes(y=log.gr,x=distir,col=drought))+
  geom_point(aes(alpha=.8))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy, labels = c("Ambient", "Reduction"))+
  facet_wrap(~year, labeller = as_labeller(labelnames.wy))+
  labs(y=" ", x="Euclidean distance to IR target", col="Precipitation 
Treatment")+
  ylim(c(-5,5))+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  guides(alpha = "none")+
  theme_ggeffects()

## ~ distance to FD traits
#wydatno21 <- wydatno21 %>% mutate(grass= ifelse(graminoid >= .51, "grassy","forby"))
distfdwy <- ggplot(wydatno21, aes(y=log.gr,x=distfd,col=drought))+
  geom_point(aes(alpha=.8))+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcolswy, labels = c("Ambient", "Reduction"))+
  facet_grid(~year, labeller = as_labeller(labelnames.wy))+
  labs(y=" ", x="Euclidean distance to FD target", col="Precipitation 
Treatment")+
  ylim(c(-5,5))+
  geom_hline(yintercept =0,col="black")+
  #stat_cor(label.y = c(c(3,3.5),c(-2.5,-2.6)))+
  #geom_hline(yintercept =0,col="black")+
  guides(alpha = "none")+
  theme_ggeffects()

#### combined figures
library(ggpubr)
wyfigtop <- ggarrange(dissboxwy,distfdwy, nrow=2,
                      common.legend = T, legend = "bottom",
                      labels = c("a","c"),label.x = .05)
wyfigbottom <-ggarrange(distdtwy,distirwy, nrow=2,
                        common.legend = T, legend = "bottom",
                        labels = c("b","d"),label.x = .05)
wyfigdrought <- ggarrange(wyfigtop,wyfigbottom, ncol=2)
wyfigdrought <- annotate_figure(wyfigdrought,
                                left="Annual growth rate")

tiff("figures/droughtfigwy.tiff", res=400, height = 5,width =8, "in",compression = "lzw")
wyfigdrought
dev.off()





################################################################

#####output of models as formatted tables for supplemental:
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
