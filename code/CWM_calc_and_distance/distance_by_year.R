#### Calculating Euclidean distances to assess dissimilarity from targets and 
#### traits minimums/maximums (because we expect lower/higher values of most traits 
#### selected to increase species qualities in a given seeding treatment).
#### Distances are calculated for all communities to every possible target (to
#### asses similarity to target + difference from other seeding trts) AND to
#### each community and its specific target (to assess ability to hit different
#### trait based targets).

# packages used
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)

# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

#### Find distance  ####
####
#### CA first
#### 
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
# irdist <- distir.ca %>% filter(trt=="ir")
# irdist <- irdist %>% 
#   unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
# # make row of targets
# irdist <- irdist %>% add_row(trt.b.y = "target", 
#                              N = quantile(traits.ca$N,.33),
#                              SRL = quantile(traits.ca$SRL,.67),
#                              RMF = quantile(traits.ca$RMF,.67)) #this should be by year tho
# irdist <- irdist %>% column_to_rownames("trt.b.y")
# #mnake into matrix
# #run dist or vegdist
# irdistmat <- vegdist(as.matrix(irdist),method = "euclidean", upper=T)#,diag=T)
# irdistmat <-as.matrix(irdistmat)
# # save only pairwise between target 
# irdistances <- as.data.frame(irdistmat["target",])
# colnames(irdistances) <- "dist"
# irdistances <- irdistances %>% rownames_to_column("trt.b.y")

## IR (min/max)
#irdist <- distir.ca %>% filter(trt=="ir")
##2021
# within each year subset the data
irdist21 <- distir.ca %>% filter(year=="2021")
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
colnames(irdistances21) <- "distir"
##2022
# within each year subset the data
irdist22 <- distir.ca %>% filter(year=="2022")
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
colnames(irdistances22) <- "distir"
##2023
# within each year subset the data
irdist23 <- distir.ca %>% filter(year=="2023")
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
colnames(irdistances23) <- "distir"
## bind dataframes
irdistances.max <- bind_rows(irdistances21,irdistances22)
irdistances.max <- bind_rows(irdistances.max,irdistances23)
irdistances.max <- irdistances.max %>% rownames_to_column("trt.b.y")

## DT
distdt.ca <- dat %>% select(c(block,trt,year,LMA,seed.mass))
distdt.ca <- merge(distdt.ca,FDdat[,c(1,2,7,10)], all.x = T)
# define trait targets
quantile(traits.ca$lma,.67) # DT
# 0.6766338 
quantile(traits.ca$seed.mass,.67) #DT
# -0.7420366
# #dtdist <- distdt.ca %>% filter(trt=="dt")
# dtdist <- distdt.ca %>% 
#   unite(trt.b.y, c(trt, block, year),sep = ".", remove=T) # make unique plot variable
# dtdist$rootdiam <- normalize(dtdist$rootdiam)
# # make row of targets
# dtdist <- dtdist %>% add_row(trt.b.y = "target", 
#                              LMA = quantile(traits.ca$LMA,.67),
#                              seed.mass = quantile(traits.ca$seed.mass,.67),
#                              rootdiam = quantile(normalize(FDdat$rootdiam),.99))
# dtdist <- dtdist %>% column_to_rownames("trt.b.y")
# #mnake into matrix
# #run dist or vegdist
# dtdistmat <- vegdist(as.matrix(dtdist),method = "euclidean", upper=T)#,diag=T)
# dtdistmat <-as.matrix(dtdistmat)
# # save only pairwise between target 
# dtdistances <- as.data.frame(dtdistmat["target",])
# colnames(dtdistances) <- "dist"
# dtdistances <- dtdistances %>% rownames_to_column("trt.b.y")


## DT (min/max)
#dtdist <- distdt.ca %>% filter(trt=="dt")
##2021
# within each year subset the data
dtdist21 <- distdt.ca %>% filter(year=="2021")
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
colnames(dtdistances21) <- "distdt"
##2022
# within each year subset the data
dtdist22 <- distdt.ca %>% filter(year=="2022")
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
colnames(dtdistances22) <- "distdt"
##2023
# within each year subset the data
dtdist23 <- distdt.ca %>% filter(year=="2023")
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
colnames(dtdistances23) <- "distdt"
## bind dataframes
dtdistances.max <- bind_rows(dtdistances21,dtdistances22)
dtdistances.max <- bind_rows(dtdistances.max,dtdistances23)
dtdistances.max <- dtdistances.max %>% rownames_to_column("trt.b.y")



## To assess similarity to highest multivariate FD we use upper 95th percentile
## of each years FD as the target to compare to. 
#FD
fddist <- FDdat[,c(1,2,9,10)] #%>% filter(trt=="fd")
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
                                 full = quantile(fddist21$full,.99)) #this should be by year tho
fddist21 <- fddist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist21),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances21 <- as.data.frame(fddistmat["target",])
colnames(fddistances21) <- "distfd"

##2022
# within each year subset the data
fddist22 <- fddist %>% filter(year=="2022")
fddist22 <- fddist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist22 <- fddist22 %>% add_row(trt.b.y = "target",
                                 full = quantile(fddist22$full,.99)) #this should be by year tho
fddist22 <- fddist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist22),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances22 <- as.data.frame(fddistmat["target",])
colnames(fddistances22) <- "distfd"

##2023
# within each year subset the data
fddist23 <- fddist %>% filter(year=="2023")
fddist23 <- fddist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist23 <- fddist23 %>% add_row(trt.b.y = "target",
                                 full = quantile(fddist23$full,.99)) #this should be by year tho
fddist23 <- fddist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist23),method = "euclidean", upper=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances23 <- as.data.frame(fddistmat["target",])
colnames(fddistances23) <- "distfd"

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
#rdist <- alldat %>% filter(trt=="rand") #%>% filter(year=="2021")

#2021
rdist21 <- alldat %>% filter(year=="2021"| year=="0")
rdist21 <- rdist21 %>% select(-c(water, Treatments,subplot))
rdist21 <- rdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
rdist21 <- rdist21 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist21),method = "euclidean", upper = T)
rdist21 <-as.matrix(randdistmat)[,-c(1:212)]
rdist21 <-rdist21[c(1:212),]
rdist21<- data.frame(
  dist=diag(as.matrix(rdist21)),
  id=colnames(rdist21))

#2022
rdist22 <- alldat %>% filter(year=="2022"| year=="0")
rdist22 <- rdist22 %>% select(-c(water, Treatments,subplot))
rdist22 <- rdist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
rdist22 <- rdist22 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist22),method = "euclidean", upper = T)
rdist22 <-as.matrix(randdistmat)[,-c(1:212)]
rdist22 <-rdist22[c(1:212),]
rdist22<- data.frame(
  dist=diag(as.matrix(rdist22)),
  id=colnames(rdist22))

#2023
rdist23 <- alldat %>% filter(year=="2023"| year=="0")
rdist23 <- rdist23 %>% select(-c(water, Treatments,subplot))
rdist23 <- rdist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
rdist23 <- rdist23 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist23),method = "euclidean", upper = T)
rdist23 <-as.matrix(randdistmat)[,-c(1:212)]
rdist23 <-rdist23[c(1:212),]
rdist23<- data.frame(
  dist=diag(as.matrix(rdist23)),
  id=colnames(rdist23))

# together 
rdistances <- bind_rows(rdist21,rdist22)
rdistances <- bind_rows(rdistances,rdist23)
rdistances <- rdistances %>% column_to_rownames("id")
colnames(rdistances) <- "distr"
rdistances <- rdistances %>% rownames_to_column("trt.b.y")
# rdistances$dist <- normalize(rdistances$dist)
# 
# irdistances$dist <- normalize(irdistances$dist)
# dtdistances$dist <- normalize(dtdistances$dist)
# 
# irdistances.max$dist <- normalize(irdistances.max$dist)
# dtdistances.max$dist <- normalize(dtdistances.max$dist)

# #### combine all dissimilarities 
# cadist <- bind_rows(dtdistances,irdistances)
# cadist <- bind_rows(cadist,fddistances)
# cadist <- bind_rows(cadist,rdistances)
# 
# #cadist$dist <- normalize(cadist$dist)
# 
# #export csv
# write.csv(cadist, "data/cwm_distances_ca.csv")

#### combine again using max in DT and IR 
cadist2 <- merge(dtdistances.max,irdistances.max)
cadist2 <- merge(cadist2,fddistances)
cadist2 <- merge(cadist2,rdistances)
# cadist2 <- bind_rows(dtdistances.max,irdistances.max)
# cadist2 <- bind_rows(cadist2,fddistances)
# cadist2 <- bind_rows(cadist2,rdistances)

## create a  column for the distance of each community to its specific target
cadist2 <- cadist2 %>% mutate(targetdist = ifelse(str_sub(trt.b.y, 1, 2)=="ir", distir,
                                                 ifelse(str_sub(trt.b.y, 1, 2)=="dt",distdt,
                                                        ifelse(str_sub(trt.b.y, 1, 2)=="fd",distfd,distr))))

#export csv
write.csv(cadist2, "data/cwm_maxdistances_ca.csv",row.names = F)



####
#### WY now
####
alldat <- read.csv("data/cwm_wy(plot).csv")
alldat$year <- as.factor(alldat$year)
alldat$block <- as.factor(alldat$block)
alldat$trt <- as.factor(alldat$trt)
FDdat <- read.csv("data/cwm_raoq_wy(plot).csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)
# CSV of species-trait combinations (for OG 25)
traits.wy <- read.csv("data/mixedgrass.csv", header=TRUE, row.names=1)
traits.wy <- traits.wy[traits.wy$use==1,] # subset use=1
traits.wy$PLSg.m2.mono <- traits.wy$PLSlbs.acre.mono * (453.59237 / 4046.8564224) #convert lb/acre to g/m2
#scale traits
traits.wy$srl = scale(log(traits.wy$srl))
traits.wy$ldmc = scale(log(traits.wy$ldmc))
traits.wy$leafn = scale(log(traits.wy$leafn))
traits.wy$lop = scale(traits.wy$lop)
traits.wy$rootdiam = scale(log(traits.wy$rootdiam))
traits.wy$sla = scale(log(traits.wy$sla))
traits.wy$rdmc = scale(log(traits.wy$rdmc))

## select relevent CWM traits per DT or IR, rows are communities
# remove predictor rows 
pred <- alldat %>% filter(year=="0")#filter preds
dat <- alldat %>% filter(year!="0")#filter preds
FDpred <- FDdat %>% filter(year=="0")#filter preds
FDdat <- FDdat %>% filter(year!="0")#filter preds

# #subset for just seeding trt of interest and calculate dissimilarity
# # first calculated for distance to target, but distance to target to be closer
# # or further from max/min of traits (ex. SRL above our target is even better
# # than hitting our target). So, distance from max/min by trt and year is also
# # now calculated below.
#
## IR
distir.wy <- dat %>% select(c(block,trt,year,leafn,srl))
distir.wy <- merge(distir.wy,FDdat[,c(2,4:6)], all.x = T)
# define trait targets
quantile(traits.wy$leafn,.25) #IR
# -0.565257
quantile(traits.wy$srl,.7557) #IR
# 0.6536232
# irdist <- distir.wy %>% filter(trt=="ir")
# irdist <- irdist %>%
#   unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
# irdist$veg <- normalize(irdist$veg)
# # make row of targets
# irdist <- irdist %>% add_row(trt.b.y = "target",
#                              leafn = quantile(traits.wy$leafn,.25),
#                              srl = quantile(traits.wy$srl,.7557),
#                              veg = quantile(normalize(FDdat$veg),.99)) #this should be by year tho
# irdist <- irdist %>% column_to_rownames("trt.b.y")
# #mnake into matrix
# #run dist or vegdist
# library(vegan)
# irdistmat <- vegdist(as.matrix(irdist),method = "euclidean", upper=T)#,diag=T)
# irdistmat <-as.matrix(irdistmat)
# # save only pairwise between target
# irdistances <- as.data.frame(irdistmat["target",])
# colnames(irdistances) <- "dist"
# irdistances <- irdistances %>% rownames_to_column("trt.b.y")

## IR (min/max)
#irdist <- distir.wy %>% filter(trt=="ir")
distir.wy$veg <- normalize(distir.wy$veg)
##2021
# within each year subset the data
irdist21 <- distir.wy %>% filter(year=="2021")
irdist21 <- irdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist21 <- irdist21 %>% add_row(trt.b.y = "target",
                                 leafn = quantile(irdist21$leafn,.05),
                                 srl = quantile(irdist21$srl,.95),
                                 veg = quantile(irdist21$veg,.95))
irdist21 <- irdist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist21),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target
irdistances21 <- as.data.frame(irdistmat["target",])
colnames(irdistances21) <- "distir"
##2022
# within each year subset the data
irdist22 <- distir.wy %>% filter(year=="2022")
irdist22 <- irdist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist22 <- irdist22 %>% add_row(trt.b.y = "target",
                                 leafn = quantile(irdist21$leafn,.05),
                                 srl = quantile(irdist21$srl,.95),
                                 veg = quantile(irdist21$veg,.95))
irdist22 <- irdist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist22),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target
irdistances22 <- as.data.frame(irdistmat["target",])
colnames(irdistances22) <- "distir"
##2023
# within each year subset the data
irdist23 <- distir.wy %>% filter(year=="2023")
irdist23 <- irdist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist23 <- irdist23 %>% add_row(trt.b.y = "target",
                                 leafn = quantile(irdist21$leafn,.05),
                                 srl = quantile(irdist21$srl,.95),
                                 veg = quantile(irdist21$veg,.95))
irdist23 <- irdist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat <- vegdist(as.matrix(irdist23),method = "euclidean", upper=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target
irdistances23 <- as.data.frame(irdistmat["target",])
colnames(irdistances23) <- "distir"
## bind dataframes
irdistances.max <- bind_rows(irdistances21,irdistances22)
irdistances.max <- bind_rows(irdistances.max,irdistances23)
irdistances.max <- irdistances.max %>% rownames_to_column("trt.b.y")


## DT
distdt.wy <- dat %>% select(c(block,trt,year,ldmc,lop))
distdt.wy <- merge(distdt.wy,FDdat[,c(1,4:6)], all.x = T)
# define trait targets
quantile(traits.wy$ldmc,.75) # DT
# 0.6766338
quantile(traits.wy$lop,.25) #DT
# -0.7420366
# dtdist <- distdt.wy %>% filter(trt=="dt")
# dtdist <- dtdist %>%
#   unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
# dtdist$rootdiam <- normalize(dtdist$rootdiam)
# # make row of targets
# dtdist <- dtdist %>% add_row(trt.b.y = "target",
#                              ldmc = quantile(traits.wy$ldmc,.75),
#                              lop = quantile(traits.wy$lop,.25),
#                              rootdiam = quantile(normalize(FDdat$rootdiam),.99)) #this should be by year tho
# dtdist <- dtdist %>% column_to_rownames("trt.b.y")
# #mnake into matrix
# #run dist or vegdist
# dtdistmat <- vegdist(as.matrix(dtdist),method = "euclidean", upper=T)#,diag=T)
# dtdistmat <-as.matrix(dtdistmat)
# # save only pairwise between target
# dtdistances <- as.data.frame(dtdistmat["target",])
# colnames(dtdistances) <- "dist"
# dtdistances <- dtdistances %>% rownames_to_column("trt.b.y")

## DT (min/max)
#dtdist <- distdt.wy %>% filter(trt=="dt")'
distdt.wy$rootdiam <- normalize(distdt.wy$rootdiam)
##2021
# within each year subset the data
dtdist21 <- distdt.wy %>% filter(year=="2021")
dtdist21 <- dtdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist21 <- dtdist21 %>% add_row(trt.b.y = "target",
                                 ldmc = quantile(dtdist21$ldmc,.95),
                                 lop = quantile(dtdist21$lop,.05),
                                 rootdiam = quantile(dtdist21$rootdiam,.95))
dtdist21 <- dtdist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist21),method = "euclidean", upper=T)#,diag=T)
dtdistmat <-as.matrix(dtdistmat)
# save only padtwise between target
dtdistances21 <- as.data.frame(dtdistmat["target",])
colnames(dtdistances21) <- "distdt"
##2022
# within each year subset the data
dtdist22 <- distdt.wy %>% filter(year=="2022")
dtdist22 <- dtdist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist22 <- dtdist22 %>% add_row(trt.b.y = "target",
                                 ldmc = quantile(dtdist21$ldmc,.95),
                                 lop = quantile(dtdist21$lop,.05),
                                 rootdiam = quantile(dtdist21$rootdiam,.95))
dtdist22 <- dtdist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist22),method = "euclidean", upper=T)#,diag=T)
dtdistmat <-as.matrix(dtdistmat)
# save only padtwise between target
dtdistances22 <- as.data.frame(dtdistmat["target",])
colnames(dtdistances22) <- "distdt"
##2023
# within each year subset the data
dtdist23 <- distdt.wy %>% filter(year=="2023")
dtdist23 <- dtdist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist23 <- dtdist23 %>% add_row(trt.b.y = "target",
                                 ldmc = quantile(dtdist21$ldmc,.95),
                                 lop = quantile(dtdist21$lop,.05),
                                 rootdiam = quantile(dtdist21$rootdiam,.95))
dtdist23 <- dtdist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist23),method = "euclidean", upper=T)
dtdistmat <-as.matrix(dtdistmat)
# save only padtwise between target
dtdistances23 <- as.data.frame(dtdistmat["target",])
colnames(dtdistances23) <- "distdt"
## bind dataframes
dtdistances.max <- bind_rows(dtdistances21,dtdistances22)
dtdistances.max <- bind_rows(dtdistances.max,dtdistances23)
dtdistances.max <- dtdistances.max %>% rownames_to_column("trt.b.y")


## To assess similarity to highest multivariate FD we use upper 95th percentile
## of each years FD as the target to compare to.
#FD
fddist <- FDdat[,c(3,4:6)]# %>% filter(trt=="fd")
#FD target (shifting annually)
quantile(normalize(fddist$full),.99)
#using max(), gives 1 in all years due to normalization, 99th quantile shows variation
subdatFD <- fddist %>% group_by(year) %>% summarize(target = quantile(normalize(full),.99))
#quantile is higher when considering less plots (i.e. only FD) and more representative of the
# comparison between FD and target FD so using subsetted data to calculate dist
# fddist <- distfd.wy %>%
#   unite(trt.b.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
fddist$full <- normalize(fddist$full)
##2021
# within each year subset the data
fddist21 <- fddist %>% filter(year=="2021")
fddist21 <- fddist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist21 <- fddist21 %>% add_row(trt.b.y = "target",
                                 full = quantile(fddist21$full,.99)) #this should be by year tho
fddist21 <- fddist21 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist21),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target
fddistances21 <- as.data.frame(fddistmat["target",])
colnames(fddistances21) <- "distfd"

##2022
# within each year subset the data
fddist22 <- fddist %>% filter(year=="2022")
fddist22 <- fddist22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist22 <- fddist22 %>% add_row(trt.b.y = "target",
                                 full = quantile(fddist22$full,.99)) #this should be by year tho
fddist22 <- fddist22 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist22),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target
fddistances22 <- as.data.frame(fddistmat["target",])
colnames(fddistances22) <- "distfd"

##2023
# within each year subset the data
fddist23 <- fddist %>% filter(year=="2023")
fddist23 <- fddist23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
fddist23 <- fddist23 %>% add_row(trt.b.y = "target",
                                 full = quantile(fddist23$full,.99)) #this should be by year tho
fddist23 <- fddist23 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist23),method = "euclidean", upper=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target
fddistances23 <- as.data.frame(fddistmat["target",])
colnames(fddistances23) <- "distfd"

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
# rdist <- alldat %>% filter(trt=="rand") %>% arrange(year,block)#%>% filter(year=="2021")
# rdist0 <- rdist %>% filter(year=="0") %>% arrange(as.numeric(block))

rdist <- alldat %>% arrange(year,block)
rdist0 <- rdist %>% filter(year=="0") %>% arrange(as.numeric(block))

#2021
rdist21 <- rdist %>% filter(year=="2021") #%>% arrange(as.numeric(block))
rdist21 <- bind_rows(rdist0,rdist21)
rdist21 <- rdist21 %>% select(-c(drought, Treatments))
rdist21 <- rdist21 %>% unite(trt.b.y, c(trt,block,year), sep = ".", remove=T) # make unique plot variable
rdist21 <- rdist21 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist21),method = "euclidean", upper = T)
rdist21 <-as.matrix(randdistmat)[,-c(1:256)]
rdist21 <-rdist21[c(1:256),]
rdist21<- data.frame(
  dist=diag(as.matrix(rdist21)),
  id=colnames(rdist21))

# rdist21 <- rdist %>% filter(year=="2021"| year=="0")
# rdist21 <- rdist21 %>% select(-c(drought, Treatments))
# rdist21 <- rdist21 %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T) # make unique plot variable
# rdist21 <- rdist21 %>% column_to_rownames("trt.b.y")
# randdistmat <- vegdist(as.matrix(rdist21),method = "euclidean", upper = T)
# rdist21 <-as.matrix(randdistmat)[,-c(1:64)]
# rdist21 <-rdist21[c(1:64),]
# rdist21.n<- data.frame(
#   dist=diag(as.matrix(rdist21)),
#   id=colnames(rdist21))
# rdist21 <- bind_rows(rdist21.n,rdist21.s)

#2022
rdist22 <- rdist %>% filter(year=="2022")
rdist22 <- bind_rows(rdist0,rdist22)
rdist22 <- rdist22 %>% select(-c(drought, Treatments))
rdist22 <- rdist22 %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T) # make unique plot variable
rdist22 <- rdist22 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist22),method = "euclidean", diag = T)
rdist22 <-as.matrix(randdistmat)[,-c(1:256)]
rdist22 <-rdist22[c(1:256),]
rdist22<- data.frame(
  dist=diag(as.matrix(rdist22)),
  id=colnames(rdist22))
# rdist22 <- rdist %>% filter(year=="2022"| year=="0") %>% filter(subplot!="s"|is.na(subplot))
# rdist22 <- rdist22 %>% select(-c(drought, Treatments))
# rdist22 <- rdist22 %>% unite(trt.b.y, c(trt, block, subplot,year), sep = ".", remove=T) # make unique plot variable
# rdist22 <- rdist22 %>% column_to_rownames("trt.b.y")
# randdistmat <- vegdist(as.matrix(rdist22),method = "euclidean", diag = T)
# rdist22 <-as.matrix(randdistmat)[,-c(1:64)]
# rdist22 <-rdist22[c(1:64),]
# rdist22.n<- data.frame(
#   dist=diag(as.matrix(rdist22)),
#   id=colnames(rdist22))
# rdist22 <- bind_rows(rdist22.n,rdist22.s)

#2023
rdist23 <- rdist %>% filter(year=="2023")
rdist23 <- bind_rows(rdist0,rdist23)
rdist23 <- rdist23 %>% select(-c(drought, Treatments))
rdist23 <- rdist23 %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T) # make unique plot variable
rdist23 <- rdist23 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist23),method = "euclidean", diag = T)
rdist23 <-as.matrix(randdistmat)[,-c(1:256)]
rdist23 <-rdist23[c(1:256),]
rdist23<- data.frame(
  dist=diag(as.matrix(rdist23)),
  id=colnames(rdist23))
# rdist23 <- rdist %>% filter(year=="2023"| year=="0") %>% filter(subplot!="s"|is.na(subplot))
# rdist23 <- rdist23 %>% select(-c(drought, Treatments))
# rdist23 <- rdist23 %>% unite(trt.b.y, c(trt, block, subplot,year), sep = ".", remove=T) # make unique plot variable
# rdist23 <- rdist23 %>% column_to_rownames("trt.b.y")
# randdistmat <- vegdist(as.matrix(rdist23),method = "euclidean", diag = T)
# rdist23 <-as.matrix(randdistmat)[,-c(1:64)]
# rdist23 <-rdist23[c(1:64),]
# rdist23.n<- data.frame(
#   dist=diag(as.matrix(rdist23)),
#   id=colnames(rdist23))
# rdist23 <- bind_rows(rdist23.n,rdist23.s)

# together
rdistances <- bind_rows(rdist21,rdist22)
rdistances <- bind_rows(rdistances,rdist23)
#rdistances <- rdistances %>% column_to_rownames("id")
colnames(rdistances) <- c("distr","trt.b.y")
# rdistances$dist <- normalize(rdistances$dist)
#
# irdistances$dist <- normalize(irdistances$dist)
# dtdistances$dist <- normalize(dtdistances$dist)
#
# irdistances.max$dist <- normalize(irdistances.max$dist)
# dtdistances.max$dist <- normalize(dtdistances.max$dist)

# #### combine all dissimilarities
# wydist <- bind_rows(dtdistances,irdistances)
# wydist <- bind_rows(wydist,fddistances)
# wydist <- bind_rows(wydist,rdistances)
#
# #wydist$dist <- normalize(wydist$dist)
#
# #export csv
# write.csv(wydist, "data/cwm_distances_wy.csv")

#### combine again using max in DT and IR
wydist2 <- merge(dtdistances.max,irdistances.max)
wydist2 <- merge(wydist2,fddistances)
wydist2 <- merge(wydist2,rdistances)
# wydist2 <- bind_rows(dtdistances.max,irdistances.max)
# wydist2 <- bind_rows(wydist2,fddistances)
# wydist2 <- bind_rows(wydist2,rdistances)

## create a  column for the distance of each community to its specific target
wydist2 <- wydist2 %>% mutate(targetdist = ifelse(str_sub(trt.b.y, 1, 2)=="ir", distir,
                                                  ifelse(str_sub(trt.b.y, 1, 2)=="dt",distdt,
                                                         ifelse(str_sub(trt.b.y, 1, 2)=="fd",distfd,distr))))

#export csv
write.csv(wydist2, "data/cwm_maxdistances_wy(plot).csv", row.names = F)





#library(patchwork)
#wyplots <- (distdt+distfd+distir+distr)
library(ggpubr)

p1 <- ggarrange(distdt,distfd,distir,distr, common.legend = T)
ggarrange(p1, alldistwy, common.legend = T, heights = c(1,2))
tiff("figures/cwm wy/alldistance_wy.tiff", res=400, height = 7,width =9, "in",compression = "lzw")
p1
dev.off()
tiff("figures/cwm ca/alldistance_ca.tiff", res=400, height = 7,width =9, "in",compression = "lzw")
p2
dev.off()

tiff("figures/cwm ca/alldistance_ca.tiff", res=400, height = 7,width =9, "in",compression = "lzw")
caplots
dev.off()

tiff("figures/cwm wy/alldistance_wy.tiff", res=400, height = 7,width =9, "in",compression = "lzw")
wyplots
dev.off()