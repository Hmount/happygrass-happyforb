#### cover analysis, are drought tolerent communities more tolerent of drought?

## packages
library(tidyverse)
library(lme4)
library(lmerTest)

## load in data and clean columns
comp.wy <- read.csv("data/hpg_total.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #remove pre-treatment data
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy.wide <- comp.wy %>% select(-c("prob","native","graminoid")) #columns to drop 
test <- comp.wy.wide %>% pivot_wider(id_cols = c("year","block","trt","subplot","drought"), 
                                names_from = "species", 
                                values_from = "cover")
comp.wy <- comp.wy %>%
  group_by(year, block, trt, subplot, drought) %>%
  filter(!any(duplicated(species) & cover == 0))
comp.wy <- comp.wy %>%
  group_by(year, block, trt, subplot, drought) %>%
  filter(!any(duplicated(species) & cover == 0.00))

comp.ca <- read.csv("data/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)]
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
#make ca long
comp.ca.long <- comp.ca %>% pivot_longer(cols = c(18:57), 
                                         names_to = "species", 
                                         values_to = "cover")


## visualize
ggplot(comp.wy.wide, aes(y=sub.tveg, x=trt, col=drought))+
  geom_boxplot()+
  facet_wrap(~year)
ggplot(comp.ca, aes(y=native.cover, x=trt, col=water))+
  geom_boxplot()+
  facet_wrap(~Year)






#model WY 
comp.wy$trt <- relevel(comp.wy$trt, ref = "r")

m1.wy <- lmer(sub.tveg ~ trt * drought + (1 | year) + (1 | block), data = comp.wy)
summary(m1.wy)
m2.wy <- lmer(sub.tveg ~ trt + drought + (1 | year), data = comp.wy)
summary(m2.wy)
anova(m1.wy,m2.wy) #better with 

var_components <- VarCorr(m1.wy)
# Assuming 'var_components' is the output from VarCorr
residual_variance <- attr(var_components, "sc")^2  # residual variance

# Total variance is the sum of residual variance and variance between groups
total_variance <- residual_variance + sum(diag(var_components))

# Calculate ICC
icc <- between_group_variance / total_variance
# Calculate ICC for each random effect
icc_random_effect1 <- var_components$year[1] / total_variance
icc_random_effect2 <- var_components$block[1] / total_variance

#repeat model for ca
comp.ca$trt <- relevel(comp.ca$trt, ref = "R")

m1.ca <- lmer(native.cover ~ trt * water + (1 | Year) + (1 | block), data = comp.ca)
summary(m1.ca)
m2.ca <- lm(native.cover ~ trt * water + (1 | Year), data = comp.ca)
summary(m2.ca)
anova(m1.ca,m2.ca) #better with 

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
