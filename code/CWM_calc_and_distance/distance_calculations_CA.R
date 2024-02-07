#### Calculating Euclidean distances to assess dissimilarity from targets and 
#### traits maximums (as high values of most traits selected should only benefit
#### species in those seeding treatments)

# packages used
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)

# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# read in CWM data, CWM_RaoQ data, and trait data
alldat <- read.csv("data/cwm_ca.csv")
alldat$year <- as.factor(alldat$year)
alldat$block <- as.factor(alldat$block)
alldat$trt <- as.factor(alldat$trt)
alldat$subplot <- as.factor(alldat$water)
FDdat <- read.csv("data/cwm_raoq_ca.csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)
FDdat$subplot <- as.factor(FDdat$water)
# CSV of species-trait combinations (for OG 25)
traits.ca <- read.csv("data/annualgrass.csv", header=TRUE, row.names=1) # CSV of species-trait combinations (for OG 25)
traits.ca <- traits.ca %>% select(-c(Asat,WUE,RLD,RTD))
# add graminoid
traits.ca$graminoid <- c(0,0,0,0,1,0,1,0,0,0, #10
                         1,1,0,1,0,1,1,1,0,0, #20
                         0,0,0,0,1,1,0,0,0,1, #30
                         1,0,0,1,1,0,1,1,0,1,1) #41
#scale traits
traits.ca$SRL = scale(log(traits.ca$SRL))
traits.ca$LMA = scale(log(traits.ca$LMA))
traits.ca$N = scale(log(traits.ca$N))
traits.ca$seed.mass = scale(traits.ca$seed.mass)
traits.ca$Rdiam = scale(log(traits.ca$Rdiam))
#traits.ca$sla = scale(log(traits.ca$sla))
traits.ca$RMF = scale(log(traits.ca$RMF))
#filter for only OG 25 species used
traits.ca <- traits.ca %>% 
  filter(Code != "AVBA") %>% 
  filter(Code != "BRNI") %>%
  filter(Code != "LUAL") %>%
  filter(Code != "MASA") %>% 
  filter(Code != "PEHE") %>% 
  filter(Code != "SACO") %>%
  filter(Code != "FEPE") %>%  # this one appears in communities, just not OG preds
  filter(Code != "FEMY") %>%  # this one appears in communities, just not OG preds
  filter(Code != "BRMA")      # this one appears in communities, just not OG preds

## select relevent CWM traits per DT or IR, rows are communities
# remove predictor rows 
pred <- alldat %>% filter(year=="0")#filter preds
dat <- alldat %>% filter(year!="0")#filter preds
FDpred <- FDdat %>% filter(year=="0")#filter preds
FDdat <- FDdat %>% filter(year!="0")#filter preds


#subset for just seeding trt of interest and calculate dissimilarity
# first calculated for distance to target, but distance to target to be closer
# or further from max/min of traits (ex. SRL above our target is even better
# than hitting our target). So, distance from max/min by trt and year is also 
# now calculated below.

## IR (targets)
distir.ca <- dat %>% select(c(block,trt,year,N,SRL,RMF))
# define trait targets
quantile(traits.ca$N,.33) #IR
# -0.565257 
quantile(traits.ca$SRL,.67) #IR
# 0.6536232 
quantile(traits.ca$RMF,.67) #IR
# 0.4271784 
irdist <- distir.ca %>% filter(trt=="ir")
irdist <- irdist %>% 
  unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
# make row of targets
irdist <- irdist %>% add_row(trt.b.y = "target", 
                             N = quantile(traits.ca$N,.33),
                             SRL = quantile(traits.ca$SRL,.67),
                             RMF = quantile(traits.ca$RMF,.67)) #this should be by year tho
irdist <- irdist %>% column_to_rownames("trt.b.y")
#mnake into matrix
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target 
irdistances <- as.data.frame(irdistmat["target",])
colnames(irdistances) <- "dist"
irdistances <- irdistances %>% rownames_to_column("trt.b.y")

## IR (min/max)
irdist <- distir.ca %>% filter(trt=="ir")
##2021
# within each year subset the data
irdist21 <- irdist %>% filter(year=="2021")
irdist21 <- irdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist21 <- irdist21 %>% add_row(trt.b.y = "target",
                                 N = quantile(irdist21$N,.05),
                                 SRL = quantile(irdist21$SRL,.95),
                                 RMF = quantile(irdist21$RMF,.95))
irdist21 <- irdist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist21),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target 
irdistances21 <- as.data.frame(irdistmat["target",])
colnames(irdistances21) <- "dist"
##2022
# within each year subset the data
irdist22 <- irdist %>% filter(year=="2022")
irdist22 <- irdist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist22 <- irdist22 %>% add_row(trt.b.y = "target",
                                 N = quantile(irdist22$N,.05),
                                 SRL = quantile(irdist22$SRL,.95),
                                 RMF = quantile(irdist22$RMF,.95))
irdist22 <- irdist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist22),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target 
irdistances22 <- as.data.frame(irdistmat["target",])
colnames(irdistances22) <- "dist"
##2023
# within each year subset the data
irdist23 <- irdist %>% filter(year=="2023")
irdist23 <- irdist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist23 <- irdist23 %>% add_row(trt.b.y = "target",
                                 N = quantile(irdist23$N,.05),
                                 SRL = quantile(irdist23$SRL,.95),
                                 RMF = quantile(irdist23$RMF,.95))
irdist23 <- irdist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist23),method = "euclidean", upper=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target 
irdistances23 <- as.data.frame(irdistmat["target",])
colnames(irdistances23) <- "dist"
## bind dataframes
irdistances.max <- bind_rows(irdistances21,irdistances22)
irdistances.max <- bind_rows(irdistances.max,irdistances23)
irdistances.max <- irdistances.max %>% rownames_to_column("trt.b.y")

## IR (bray-curtis)
irdist <- distir.ca %>% filter(trt=="ir")
irdist <- irdist %>% 
  unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
# make row of targets
irdist <- irdist %>% add_row(trt.b.y = "target", 
                             N = quantile(traits.ca$N,.33),
                             SRL = quantile(traits.ca$SRL,.67),
                             RMF = quantile(traits.ca$RMF,.67)) #this should be by year tho
irdist <- irdist %>% column_to_rownames("trt.b.y")
#mnake into matrix
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target 
irdistances.bc <- as.data.frame(irdistmat["target",])
colnames(irdistances.bc) <- "dist"
irdistances.bc <- irdistances.bc %>% rownames_to_column("trt.b.y")

## DT
distdt.ca <- dat %>% select(c(block,trt,year,LMA,seed.mass))
distdt.ca <- merge(distdt.ca,FDdat[,c(1,2,7,10)], all.x = T)
# define trait targets
quantile(traits.ca$lma,.67) # DT
# 0.6766338 
quantile(traits.ca$seed.mass,.67) #DT
# -0.7420366
dtdist <- distdt.ca %>% filter(trt=="dt")
dtdist <- dtdist %>% 
  unite(trt.b.y, c(trt, block, year),sep = ".", remove=T) # make unique plot variable
dtdist$rootdiam <- normalize(dtdist$rootdiam)
# make row of targets
dtdist <- dtdist %>% add_row(trt.b.y = "target", 
                             LMA = quantile(traits.ca$LMA,.67),
                             seed.mass = quantile(traits.ca$seed.mass,.67),
                             rootdiam = quantile(normalize(FDdat$rootdiam),.99))
dtdist <- dtdist %>% column_to_rownames("trt.b.y")
#mnake into matrix
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist),method = "euclidean", upper=T)#,diag=T)
dtdistmat <-as.matrix(dtdistmat)
# save only pairwise between target 
dtdistances <- as.data.frame(dtdistmat["target",])
colnames(dtdistances) <- "dist"
dtdistances <- dtdistances %>% rownames_to_column("trt.b.y")


## DT (min/max)
dtdist <- distdt.ca %>% filter(trt=="dt")
##2021
# within each year subset the data
dtdist21 <- dtdist %>% filter(year=="2021")
dtdist21 <- dtdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist21 <- dtdist21 %>% add_row(trt.b.y = "target", 
                                 LMA = quantile(dtdist21$LMA,.95),
                                 seed.mass = quantile(dtdist21$seed.mass,.95),
                                 rootdiam = quantile(normalize(dtdist21$rootdiam),.95))
dtdist21 <- dtdist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist21),method = "euclidean", upper=T)#,diag=T)
dtdistmat <-as.matrix(dtdistmat)
# save only padtwise between target 
dtdistances21 <- as.data.frame(dtdistmat["target",])
colnames(dtdistances21) <- "dist"
##2022
# within each year subset the data
dtdist22 <- dtdist %>% filter(year=="2022")
dtdist22 <- dtdist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist22 <- dtdist22 %>% add_row(trt.b.y = "target", 
                                 LMA = quantile(dtdist22$LMA,.95),
                                 seed.mass = quantile(dtdist22$seed.mass,.95),
                                 rootdiam = quantile(normalize(dtdist22$rootdiam),.95))
dtdist22 <- dtdist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist22),method = "euclidean", upper=T)#,diag=T)
dtdistmat <-as.matrix(dtdistmat)
# save only padtwise between target 
dtdistances22 <- as.data.frame(dtdistmat["target",])
colnames(dtdistances22) <- "dist"
##2023
# within each year subset the data
dtdist23 <- dtdist %>% filter(year=="2023")
dtdist23 <- dtdist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist23 <- dtdist23 %>% add_row(trt.b.y = "target", 
                                 LMA = quantile(dtdist23$LMA,.95),
                                 seed.mass = quantile(dtdist23$seed.mass,.95),
                                 rootdiam = quantile(normalize(dtdist23$rootdiam),.95))
dtdist23 <- dtdist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist23),method = "euclidean", upper=T)
dtdistmat <-as.matrix(dtdistmat)
# save only padtwise between target 
dtdistances23 <- as.data.frame(dtdistmat["target",])
colnames(dtdistances23) <- "dist"
## bind dataframes
dtdistances.max <- bind_rows(dtdistances21,dtdistances22)
dtdistances.max <- bind_rows(dtdistances.max,dtdistances23)
dtdistances.max <- dtdistances.max %>% rownames_to_column("trt.b.y")



## To assess similarity to highest multivariate FD we use upper 95th percentile
## of each years FD as the target to compare to. 
#FD
fddist <- FDdat[,c(1,2,9,10)] %>% filter(trt=="fd")
#FD target (shifting annually)
quantile(normalize(fddist$full),.95) 
#using max(), gives 1 in all years due to normalization, 99th quantile shows variation
subdatFD <- fddist %>% group_by(year) %>% summarize(target = quantile(normalize(full),.95))
#quantile is higher when considering less plots (i.e. only FD) and more representative of the 
# comparison between FD and target FD so using subsetted data to calculate dist
# fddist <- distfd.ca %>% 
#   unite(trt.b.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
fddist$full <- normalize(fddist$full)
##2021
# within each year subset the data
fddist21 <- fddist %>% filter(year=="2021")
fddist21 <- fddist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist21 <- fddist21 %>% add_row(trt.b.y = "target",
                                 full = quantile(normalize(fddist21$full),.99)) #this should be by year tho
fddist21 <- fddist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist21),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances21 <- as.data.frame(fddistmat["target",])
colnames(fddistances21) <- "dist"

##2022
# within each year subset the data
fddist22 <- fddist %>% filter(year=="2022")
fddist22 <- fddist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist22 <- fddist22 %>% add_row(trt.b.y = "target",
                                 full = quantile(normalize(fddist22$full),.99)) #this should be by year tho
fddist22 <- fddist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist22),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances22 <- as.data.frame(fddistmat["target",])
colnames(fddistances22) <- "dist"

##2023
# within each year subset the data
fddist23 <- fddist %>% filter(year=="2023")
fddist23 <- fddist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist23 <- fddist23 %>% add_row(trt.b.y = "target",
                                 full = quantile(normalize(fddist23$full),.99)) #this should be by year tho
fddist23 <- fddist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist23),method = "euclidean", upper=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances23 <- as.data.frame(fddistmat["target",])
colnames(fddistances23) <- "dist"

#bind dataframes
fddistances <- bind_rows(fddistances21,fddistances22)
fddistances <- bind_rows(fddistances,fddistances23)
fddistances <- fddistances %>% rownames_to_column("trt.b.y")


## To assess similarity to random so it can be included in models, we consider the distances 
## in CWM from target seeding CWM (could also look at community abundance/ composition from 
## our expected probabilities from seeding, but the cwm matches other dist measures.)
## this essentially uses the random community as a measure of how well be were at making any 
## CWM we sought out to make by seeding and allows us to continue using random as a control. 
# Rand

# within each year subset the data
rdist <- alldat %>% filter(trt=="rand") #%>% filter(year=="2021")

#2021
rdist21 <- rdist %>% filter(year=="2021"| year=="0")
rdist21 <- rdist21 %>% select(-c(water, Treatments,subplot))
rdist21 <- rdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
rdist21 <- rdist21 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist21),method = "euclidean", upper = T)
rdist21 <-as.matrix(randdistmat)[,-c(1:53)]
rdist21 <-rdist21[c(1:53),]
rdist21<- data.frame(
  dist=diag(as.matrix(rdist21)),
  id=colnames(rdist21))

#2022
rdist22 <- rdist %>% filter(year=="2022"| year=="0")
rdist22 <- rdist22 %>% select(-c(water, Treatments,subplot))
rdist22 <- rdist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
rdist22 <- rdist22 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist22),method = "euclidean", upper = T)
rdist22 <-as.matrix(randdistmat)[,-c(1:53)]
rdist22 <-rdist22[c(1:53),]
rdist22<- data.frame(
  dist=diag(as.matrix(rdist22)),
  id=colnames(rdist22))

#2023
rdist23 <- rdist %>% filter(year=="2023"| year=="0")
rdist23 <- rdist23 %>% select(-c(water, Treatments,subplot))
rdist23 <- rdist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
rdist23 <- rdist23 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist23),method = "euclidean", upper = T)
rdist23 <-as.matrix(randdistmat)[,-c(1:53)]
rdist23 <-rdist23[c(1:53),]
rdist23<- data.frame(
  dist=diag(as.matrix(rdist23)),
  id=colnames(rdist23))

# together 
rdistances <- bind_rows(rdist21,rdist22)
rdistances <- bind_rows(rdistances,rdist23)
rdistances <- rdistances %>% column_to_rownames("id")
colnames(rdistances) <- "dist"
rdistances <- rdistances %>% rownames_to_column("trt.b.y")
# rdistances$dist <- normalize(rdistances$dist)
# 
# irdistances$dist <- normalize(irdistances$dist)
# dtdistances$dist <- normalize(dtdistances$dist)
# 
# irdistances.max$dist <- normalize(irdistances.max$dist)
# dtdistances.max$dist <- normalize(dtdistances.max$dist)

#### combine all dissimilarities 
cadist <- bind_rows(dtdistances,irdistances)
cadist <- bind_rows(cadist,fddistances)
cadist <- bind_rows(cadist,rdistances)

#cadist$dist <- normalize(cadist$dist)

#export csv
write.csv(cadist, "data/cwm_distances_ca.csv")

#### combine again using max in DT and IR 
cadist2 <- bind_rows(dtdistances.max,irdistances.max)
cadist2 <- bind_rows(cadist2,fddistances)
cadist2 <- bind_rows(cadist2,rdistances)

#export csv
write.csv(cadist2, "data/cwm_maxdistances_ca.csv")



#### Figures
#combine with CWM to plot
cwm.ca <- read.csv("data/cwm_ca.csv")# California CWM data
cwm.ca <- cwm.ca %>% filter(year != "0") #remove predicted/target communities
cwm.ca$year <- as.factor(cwm.ca$year)
cwm.ca <- cwm.ca %>% mutate(yrorder = ifelse(year=="2021","1",
                                      ifelse(year=="2022","2",
                                      ifelse(year=="2023","3","0")))) #make new sequence column
cwm.ca$yrorder <- as.numeric(cwm.ca$yrorder)
cwm.ca <- cwm.ca %>% 
  mutate(plot = paste(block, trt, sep = ".")) #add plot ID column (but give NA to target/predicted communities)
cwm.ca <- cwm.ca %>% filter(year != "0") #remove predicted/target communities
# also combine CWM_distances dataframe to master df 
#break apart distances ID to make wider and merge together
cadist2 <- separate(cadist2, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
cadist2 <- cadist2 %>% filter(trt!="target")
#cadist <- cadist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allca <- merge(cwm.ca,cadist2, by=c("year","trt","block"), all=T)
allca$trt <- as.factor(allca$trt)
allca$water <- as.factor(allca$water)

#make dought column
allca <- allca %>% mutate(drought = ifelse(water=="0.5","drt","cntl"))

#set reference levels for modelling
allca$trt <- relevel(allca$trt, ref = "rand") #make random communities the reference level
allca$drought <- as.factor(allca$drought)
allca$drought <- relevel(allca$drought, ref = "cntl") #make random communities the reference level

ggplot(allca,(aes(x=dist)))+
  geom_histogram()+
  facet_wrap(~trt) #looks fiarly normal within groups
hist(allca$dist)
shapiro.test(allca$dist) #below .05 = not normal!
# hist(sqrt(allca$dist))
# shapiro.test(sqrt(allca$dist)) #still not normal, but closer
# allca$dist_tran <- sqrt(allca$dist)
### for now keeping dist untransformed 

#transfomed
# summary(t <- aov(dist_tran~trt*drought*year, allca))
# tuktest <- TukeyHSD(t)
# multcompView::multcompLetters4(t,tuktest)
# ggplot(allca, aes(y=dist_tran, x=trt, fill=drought))+
#   geom_boxplot()+
#   #geom_smooth(method="lm")+
#   facet_wrap(~year, scales="fixed")+
#   labs(x=" ",y="Distance from target (normalized)")+ #, fill="drought treatment")+
#   theme_classic()+
#   ylim(0,3)

#w/o transform
summary(t <- aov(dist~trt*drought*year, allca))
tuktest <- TukeyHSD(t)
multcompView::multcompLetters4(t,tuktest)
ggplot(allca, aes(y=dist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Distance from target (normalized)")+ #, fill="drought treatment")+
  theme_classic()+
  ylim(0,4)

letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought:year'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\2", rownames(letterstest)))
letterstest$year <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\3", rownames(letterstest)))
test <- allca %>% group_by(year, drought, trt) %>% summarise(yposition = quantile(dist,.8))
test <- merge(letterstest,test, by = c("year", "drought", "trt"))
test2 <- merge(test,allca, by = c("year", "drought", "trt"), all=T)

droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

#export for short report
tiff("figures/cwm ca/distanceplot_ca.tiff", res=400, height = 4,width =8.5, "in",compression = "lzw")
ggplot(test2, aes(y=dist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  scale_fill_manual(values = droughtcolsca)+
  #facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Distance from target")+ #, fill="drought treatment")+
  theme_classic()+
  ylim(0,4)
dev.off()

#combined with distance by year on other r script

letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- allca %>% group_by(drought, trt) %>% summarise(yposition = quantile(dist,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,allca, by = c("drought", "trt"), all=T)

droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color


alldistca <- ggplot(test2, aes(y=dist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  scale_fill_manual(values = droughtcolswy)+
  #facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Distance from target")+ #, fill="drought treatment")+
  theme_classic()+
  ylim(0,4)
