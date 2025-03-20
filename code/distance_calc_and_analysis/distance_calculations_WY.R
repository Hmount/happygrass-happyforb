#### Calculating Euclidean distances to assess dissimilarity of community 
#### CWM traits and Roa/Q traits from trait targets in each plot.
#### Traits are still considered as achieving the target if they exceed the 
#### upper quantile or are below the minimum quantile (because minimum + 
#### maximum trait values should achieve their target even better than 
#### upper/lower quantile, respectively.). Targets are calculated within 
#### each year, and the functionally diverse goal was allowed to vary by 
#### year because different maximum FD was possible at each site each year 
#### by virtue of the species pool present. 
#### Wyoming ####

# packages used
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)

# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

#### load and clean data
alldat <- read.csv("data/cwm_wy(plot).csv")
alldat$year <- as.factor(alldat$year)
alldat$block <- as.factor(alldat$block)
alldat$trt <- as.factor(alldat$trt)
FDdat <- read.csv("data/cwm_raoq_wy(plot).csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)
# CSV of species-trait combinations (for OG 25)
traits.wy <- read.csv("data/trait_data/mixedgrass.csv", header=TRUE, row.names=1)
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

#### Find distance
#### subset for just seeding trt of interest and calculate dissimilarity
#### First distance to exact quantile-based target was used (commented out),
#### but instead we now find distance to the max/min of trait quantiles 
#### because communities still achieved if the target is above the upper quantile
#### or below the lower qualtile (ex. SRL above our target is even better
#### than hitting our target). We use the max/min distances for subsequent analyses.

## IR
distir.wy <- dat %>% select(c(block,trt,year,leafn,srl))
distir.wy <- merge(distir.wy,FDdat[,c(2,4:6)], all.x = T)
# define trait targets
quantile(traits.wy$leafn,.25) #IR
# -0.565257
quantile(traits.wy$srl,.7557) #IR
# 0.6536232

## (min/max)
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


#### for RC:
#### We cannot compare the random community to a trait-based target because we
#### established these communities with no target (using a log normal distribution).
#### Instead, we use the CWM traits that the random control community should have 
#### produced based on the relative abundance of species in the seed mix and 
#### compare that to the realized CWM's for all six traits each year. 
#### This is useful to calculate for RC because it provides a control relationship between
#### functional change and taxonomic change (bray-curtis) like we examine with the other 
#### treatment groups. 

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

# together
rdistances <- bind_rows(rdist21,rdist22)
rdistances <- bind_rows(rdistances,rdist23)
colnames(rdistances) <- c("distr","trt.b.y")

#### combine using min/max in DT and IR
wydist2 <- merge(dtdistances.max,irdistances.max)
wydist2 <- merge(wydist2,fddistances)
wydist2 <- merge(wydist2,rdistances)

## create a  column for the distance of each community to its specific target
wydist2 <- wydist2 %>% mutate(targetdist = ifelse(str_sub(trt.b.y, 1, 2)=="ir", distir,
                                                  ifelse(str_sub(trt.b.y, 1, 2)=="dt",distdt,
                                                         ifelse(str_sub(trt.b.y, 1, 2)=="fd",distfd,distr))))

#export csv
write.csv(wydist2, "data/cwm_maxdistances_wy(plot).csv", row.names = F)
