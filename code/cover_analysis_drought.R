#### cover analysis, are drought tolerent communities more tolerent of drought?

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## load in data, clean and modify columns
# WY
comp.wy <- read.csv("data/hpg_total.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
# add invaded locations as 0/1 variable (2023 only)
wy.invaded.23 <- read.csv("data/invasion_loc23.csv") # load data
wy.invaded.23$invaded <-tolower(wy.invaded.23$invaded)
comp.wy <- comp.wy %>%
  mutate(invaded = ifelse(year == 2023, as.integer(paste(block, trt, subplot) %in% paste(wy.invaded.23$block, wy.invaded.23$trt, wy.invaded.23$invaded)), NA_integer_))
rm(wy.invaded.23)# remove invasion data
# make unique plot variable
comp.wy <- comp.wy %>% unite(plot, c(block, trt, subplot), sep = ".", remove=F) 
# calculate cover native per subplot and add to data
fornativecover <- comp.wy %>% filter(species!="BG"&
                        species!="Litter"&
                        native == "N") %>% #only native live veg
  group_by(year,block,trt,subplot) %>% 
  summarize(nativecov = sum(cover, na.rm=T)) #summarize total live veg per subplot
comp.wy <- merge(comp.wy,fornativecover, all.x = T)
# calculate cover native per PLOT (averaged per subplot)
forplotcover <- comp.wy %>% filter(species!="BG"&
                                       species!="Litter"&
                                       native == "N") %>% #only native live veg
  group_by(year,block,trt) %>% 
  summarize(nativecov.plotmean = sum(cover, na.rm=T)/2) #summarize total live veg per PLOT
comp.wy <- merge(comp.wy,fornativecover, all.x = T)
# make wide for analysis and matching CA data
comp.wy.wide <- comp.wy %>% select(-c("prob","native","graminoid")) #columns to drop 
comp.wy.wide <- comp.wy.wide %>% pivot_wider(id_cols = c("year","block","trt","subplot","drought",
                                                         "nativecov","BG", "Litter","plot", "invaded"), 
                                names_from = "species", 
                                values_from = "cover")
comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data

# CA
comp.ca <- read.csv("data/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)] #remove empty columns
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
#make ca long
#comp.ca.long <- comp.ca %>% pivot_longer(cols = c(18:57), 
#                                         names_to = "species", 
#                                         values_to = "cover")


### data summary
# look at response variable in each dataset
hist(comp.wy.wide$nativecov) # native cover in WY is left skewed
hist(sqrt(comp.wy.wide$nativecov)) # good transformation for normalizing data

hist(comp.ca$native.cover) # native cover in CA is not skewed
hist(sqrt(comp.ca$native.cover)) # but transform to match


### visualize
# WY
droughtcolswy <- c("0"="skyblue", "1"="red") #create variable for color
ggplot(comp.wy.wide, aes(y=nativecov, x=trt, fill=drought))+
  geom_boxplot()+
  scale_fill_manual(values = droughtcolswy)+
  facet_wrap(~year)

# CA
#droughtcolsca <- c(".5"="skyblue", "1.25"="red")
ggplot(comp.ca, aes(y=native.cover, x=trt, fill=water))+
  geom_boxplot()+
  #scale_fill_manual(values = droughtcolsca)+
  facet_wrap(~Year)


### modelling native community cover as a function of drought and seeding community
# WY 
comp.wy.wide$trt <- relevel(comp.wy.wide$trt, ref = "r") #make random communities the reference level

m0.wy <- lmer(sqrt(nativecov) ~ trt * drought + (1 | year) + (1 | block), data = comp.wy.wide)
summary(m0.wy) # simplest model, but doesn't account for invasion in 2023
m1.wy <- lmer(sqrt(nativecov) ~ trt * drought + (1 | year) + (1 | block) + (1|plot), data = comp.wy.wide)
summary(m1.wy) # plot can account for plot effect
anova(m0.wy,m1.wy) # plot does not improve model statistically, but should be retained

subun <- comp.wy.wide %>% filter(invaded == 0 | is.na(invaded)) #subset data for uninvaded only
m2.wy <- lmer(sqrt(nativecov) ~ trt * drought  + (1 | year) + (1 | block), data = subun) #without invaded
summary(m2.wy) #fd and fd*drought drop out

## Best model:
summary(m1.wy)
anova(m1.wy)

# assess random effects
var_components <- VarCorr(m1.wy)
# Assuming 'var_components' is the output from VarCorr
residual_variance <- attr(var_components, "sc")^2  # residual variance
# Total variance is the sum of residual variance and variance between groups
total_variance <- residual_variance + sum(diag(var_components))
# Calculate ICC
icc <- between_group_variance / total_variance
# Calculate ICC for each random effect
icc_random_effect1 <- var_components$year[1] / total_variance #.582
icc_random_effect2 <- var_components$block[1] / total_variance #.066


# CA
comp.ca$trt <- relevel(comp.ca$trt, ref = "R") # random as refernce level
comp.ca$water <- relevel(comp.ca$water, ref = "1.25") # water as refernece level (drought = treatment)

m0.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year) + (1 | block) + (1| plot), data = comp.ca)
summary(m0.ca) # this matches WY model, but not converging (can't have block + plot)
m1.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year) + (1 | plot), data = comp.ca)
summary(m1.ca) 
# plot probably makes more sense given design and more variance is explained in plot than block
# but do these models match? WHY?

# assess random effects
var_componentsca <- VarCorr(m1.ca)
# Assuming 'var_components' is the output from VarCorr
residual_varianceca <- attr(var_componentsca, "sc")^2  # residual variance
# Total variance is the sum of residual variance and variance between groups
total_varianceca <- residual_varianceca + sum(diag(var_componentsca))
# Calculate ICC
icc <- between_group_varianceca / total_varianceca
# Calculate ICC for each random effect
icc_random_effect1ca <- var_componentsca$Year[1] / total_varianceca
icc_random_effect2ca <- var_componentsca$block[1] / total_varianceca
