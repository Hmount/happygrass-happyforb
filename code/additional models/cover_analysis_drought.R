#### cover analysis, are drought tolerent communities more tolerent of drought?

## packages
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

## load in data, clean and modify columns
# WY
comp.wy <- read.csv("data/comp_wy.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy$nativecov <- comp.wy$nativecov/100  # make native live veg % a proportion to match CA data
comp.wy$plot.tveg <- comp.wy$plot.tveg/100  # make native live veg % a proportion to match CA data

# CA
comp.ca <- read.csv("data/comp_ca.csv") #read in California comp data
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)


### data summary
# look at response variable in each dataset
# WY
hist(comp.wy$nativecov) # native cover in WY is left skewed
hist(sqrt(comp.wy$nativecov)) # good transformation for normalizing data
hist(comp.wy$plot.tveg) # native cover (mean plot) in WY is also skewed
hist(sqrt(comp.wy$plot.tveg)) # same transformation works for normalizing data
# CA
hist(comp.ca$native.cover) # native cover in CA is not skewed
hist(sqrt(comp.ca$native.cover)) # but transform to match


### visualize
# WY
droughtcolswy <- c("0"="skyblue", "1"="red") #create variable for color
ggplot(comp.wy.wide, aes(y=nativecov, x=trt, fill=drought))+ #subplot 
  geom_boxplot()+
  scale_fill_manual(values = droughtcolswy)+
  facet_wrap(~year)
ggplot(allwy, aes(y=nativecov, x=trt, fill=drought))+ #plot
  geom_boxplot()+
  scale_fill_manual(values = droughtcols)+
  facet_wrap(~year)

# CA
#droughtcolsca <- c(".5"="skyblue", "1.25"="red")
ggplot(comp.ca, aes(y=native.cover, x=trt, fill=water))+
  geom_boxplot()+
  #scale_fill_manual(values = droughtcolsca)+
  facet_wrap(~Year)


### modelling native community cover as a function of drought and seeding community
#### WY ####
comp.wy.wide$trt <- relevel(comp.wy.wide$trt, ref = "r") #make random communities the reference level
# begin modelling with subplots, however these models are ultimately dropped in favor of 
# models of cover at the plot level that better match the CA data.
m0.wy <- lmer(sqrt(nativecov) ~ trt * drought + (1 | year) + (1 | block), data = comp.wy.wide)
summary(m0.wy) # simplest model, but doesn't account for invasion in 2023
m1.wy <- lmer(sqrt(nativecov) ~ trt * drought + (1 | year) + (1 | block) + (1|plot), data = comp.wy.wide)
summary(m1.wy) # plot can account for plot effect
anova(m0.wy,m1.wy) # plot does not improve model statistically, but should be retained

subun <- comp.wy.wide %>% filter(invaded == 0 | is.na(invaded)) #subset data for uninvaded only
m2.wy <- lmer(sqrt(nativecov) ~ trt * drought  + (1 | year) + (1 | block), data = subun) #without invaded
summary(m2.wy) #fd and fd*drought drop out

# Best model (at subplot level):
summary(m1.wy)
anova(m1.wy)

## test for differences when using plot instead of subplot (more like CA)
m1.wy.plot <- lmer(sqrt(plot.tveg) ~ trt * drought + (1 | year) + (1|plot), data = comp.wy)
summary(m1.wy.plot) # plot or block, but not both (plot matches subplot model)
# # appears to be little to no difference when using plot in WY 
# # plus, using plot matches what I was using for smaller sample size CA models
# anova(m1.wy.plot,m1.wy.plotb)
# AIC(m1.wy.plot)

# new best model (plot level)
summary(m1.wy.plot)
anova(m1.wy.plot)
emm.wy <- emmeans(m1.wy.plot, c("trt","drought"))
pairs(emm.wy) # significant differences between most groups (but not much in drought)

# assess random effects
var_components <- VarCorr(m1.wy.plot)
# Assuming 'var_components' is the output from VarCorr
residual_variance <- attr(var_components, "sc")^2  # residual variance
# Total variance is the sum of residual variance and variance between groups
total_variance <- residual_variance + sum(diag(var_components))
# Calculate ICC
icc <- between_group_variance / total_variance
# Calculate ICC for each random effect
icc_random_effect1 <- var_components$year[1] / total_variance #.737
icc_random_effect2 <- var_components$plot[1] / total_variance #.056

## What about adding CWM to the model?
m2.wy.plot <- lmer(sqrt(plot.tveg) ~ trt * drought * CWM + (1 | year) + (1|plot), data = comp.wy)


#### CA ####
comp.ca$trt <- relevel(comp.ca$trt, ref = "R") # random as reference level
comp.ca$water <- relevel(comp.ca$water, ref = "1.25") # water as reference level (drought = treatment)

m0.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year) + (1 | block) + (1| plot), data = comp.ca)
summary(m0.ca) # this matches WY model, but not converging (can't have block + plot)
m1.ca <- lmer(sqrt(native.cover) ~ trt * water + (1 | Year) + (1 | plot), data = comp.ca)
summary(m1.ca) 
# # plot probably makes more sense given design and more variance is explained in plot than block
# anova(m1.ca,m1.cab)

# best model
summary(m1.ca) 
anova(m1.ca)
emm.ca <- emmeans(m1.ca, c("trt","water"))
pairs(emm.ca) # not really significant differences between any groups

# assess random effects
var_componentsca <- VarCorr(m1.ca)
# Assuming 'var_components' is the output from VarCorr
residual_varianceca <- attr(var_componentsca, "sc")^2  # residual variance
# Total variance is the sum of residual variance and variance between groups
total_varianceca <- residual_varianceca + sum(diag(var_componentsca))
# # Calculate ICC
# icc <- between_group_varianceca / total_varianceca
# Calculate ICC for each random effect
icc_random_effect1ca <- var_componentsca$Year[1] / total_varianceca # .035
icc_random_effect2ca <- var_componentsca$plot[1] / total_varianceca # .268
