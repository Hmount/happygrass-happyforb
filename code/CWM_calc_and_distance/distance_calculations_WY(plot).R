#### Calculating Euclidean distances to assess dissimilarity from targets and 
#### traits maximums (as high values of most traits selected should only benefit
#### species in those seeding treatments)

########## This is a repeat of "distance_calculations_WY.R" at PLOT level ############

# packages used
library(tidyverse)
library(lme4)
library(lmerTest)
library(vegan)

# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

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
irdist <- distir.wy %>% filter(trt=="ir")
irdist <- irdist %>%
  unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
irdist$veg <- normalize(irdist$veg)
# make row of targets
irdist <- irdist %>% add_row(trt.b.y = "target",
                             leafn = quantile(traits.wy$leafn,.25),
                             srl = quantile(traits.wy$srl,.7557),
                             veg = quantile(normalize(FDdat$veg),.99)) #this should be by year tho
irdist <- irdist %>% column_to_rownames("trt.b.y")
#mnake into matrix
#run dist or vegdist
library(vegan)
irdistmat <- vegdist(as.matrix(irdist),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target
irdistances <- as.data.frame(irdistmat["target",])
colnames(irdistances) <- "dist"
irdistances <- irdistances %>% rownames_to_column("trt.b.y")

## IR (min/max)
irdist <- distir.wy %>% filter(trt=="ir")
##2021
# within each year subset the data
irdist21 <- irdist %>% filter(year=="2021")
irdist21 <- irdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
irdist21 <- irdist21 %>% add_row(trt.b.y = "target",
                                 leafn = quantile(irdist21$leafn,.05),
                                 srl = quantile(irdist21$srl,.95),
                                 veg = quantile(normalize(irdist21$veg),.95))
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
                                 leafn = quantile(irdist21$leafn,.05),
                                 srl = quantile(irdist21$srl,.95),
                                 veg = quantile(normalize(irdist21$veg),.95))
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
                                 leafn = quantile(irdist21$leafn,.05),
                                 srl = quantile(irdist21$srl,.95),
                                 veg = quantile(normalize(irdist21$veg),.95))
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


## DT
distdt.wy <- dat %>% select(c(block,trt,year,ldmc,lop))
distdt.wy <- merge(distdt.wy,FDdat[,c(1,4:6)], all.x = T)
# define trait targets
quantile(traits.wy$ldmc,.75) # DT
# 0.6766338
quantile(traits.wy$lop,.25) #DT
# -0.7420366
dtdist <- distdt.wy %>% filter(trt=="dt")
dtdist <- dtdist %>%
  unite(trt.b.y, c(trt, block, year), sep = ".", remove=T) # make unique plot variable
dtdist$rootdiam <- normalize(dtdist$rootdiam)
# make row of targets
dtdist <- dtdist %>% add_row(trt.b.y = "target",
                             ldmc = quantile(traits.wy$ldmc,.75),
                             lop = quantile(traits.wy$lop,.25),
                             rootdiam = quantile(normalize(FDdat$rootdiam),.99)) #this should be by year tho
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
dtdist <- distdt.wy %>% filter(trt=="dt")
##2021
# within each year subset the data
dtdist21 <- dtdist %>% filter(year=="2021")
dtdist21 <- dtdist21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# attach row of targets
dtdist21 <- dtdist21 %>% add_row(trt.b.y = "target",
                                 ldmc = quantile(dtdist21$ldmc,.95),
                                 lop = quantile(dtdist21$lop,.05),
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
                                 ldmc = quantile(dtdist21$ldmc,.95),
                                 lop = quantile(dtdist21$lop,.05),
                                 rootdiam = quantile(normalize(dtdist21$rootdiam),.95))
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
                                 ldmc = quantile(dtdist21$ldmc,.95),
                                 lop = quantile(dtdist21$lop,.05),
                                 rootdiam = quantile(normalize(dtdist21$rootdiam),.95))
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
fddist <- FDdat[,c(3,4:6)] %>% filter(trt=="fd")
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
rdist <- alldat %>% filter(trt=="rand") %>% arrange(year,block)#%>% filter(year=="2021")
rdist0 <- rdist %>% filter(year=="0") %>% arrange(as.numeric(block))

#2021
rdist21 <- rdist %>% filter(year=="2021") #%>% arrange(as.numeric(block))
rdist21 <- bind_rows(rdist0,rdist21)
rdist21 <- rdist21 %>% select(-c(drought, Treatments))
rdist21 <- rdist21 %>% unite(trt.b.y, c(trt,block,year), sep = ".", remove=T) # make unique plot variable
rdist21 <- rdist21 %>% column_to_rownames("trt.b.y")
randdistmat <- vegdist(as.matrix(rdist21),method = "euclidean", upper = T)
rdist21 <-as.matrix(randdistmat)[,-c(1:64)]
rdist21 <-rdist21[c(1:64),]
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
rdist22 <-as.matrix(randdistmat)[,-c(1:64)]
rdist22 <-rdist22[c(1:64),]
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
rdist23 <-as.matrix(randdistmat)[,-c(1:64)]
rdist23 <-rdist23[c(1:64),]
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
colnames(rdistances) <- c("dist","trt.b.y")
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
wydist2 <- bind_rows(dtdistances.max,irdistances.max)
wydist2 <- bind_rows(wydist2,fddistances)
wydist2 <- bind_rows(wydist2,rdistances)

#export csv
write.csv(wydist2, "data/cwm_maxdistances_wy(plot).csv")



#### Figures

#make new sequence column
#add plot ID column (but give NA to target/predicted communities)

# # combine to master df (remove spp columns for now)
# allwy <- merge(comp.wy[,-c(11:66)],cwm.wy, by=c("year","trt","block","subplot"))
# allwy$trt <- as.factor(allwy$trt)


#combine with CWM to plot
cwm.wy <- read.csv("data/cwm_wy.csv")# Wyoming CWM data
cwm.wy$year <- as.factor(cwm.wy$year)
cwm.wy <- cwm.wy %>% mutate(yrorder = ifelse(year=="2021","1",
                                             ifelse(year=="2022","2",
                                                    ifelse(year=="2023","3","0")))) #make new sequence column
cwm.wy$yrorder <- as.numeric(cwm.wy$yrorder)
cwm.wy <- cwm.wy %>% 
  mutate(plot = paste(block, trt, year, sep = ".")) #add plot ID column (but give NA to target/predicted communities)
cwm.wy <- cwm.wy %>% filter(year != "0") #remove predicted/target communities
# also combine CWM_distances dataframe to master df 
#break apart distances ID to make wider and merge together
wydist2 <- separate(wydist2, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
wydist2 <- wydist2 %>% filter(trt!="target")
#wydist <- wydist %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
allwy <- merge(cwm.wy,wydist2, by=c("year","trt","block"), all.x=T)
allwy$trt <- as.factor(allwy$trt)
allwy$drought <- as.factor(allwy$drought)

#set reference levels for modelling
allwy$trt <- relevel(allwy$trt, ref = "rand") #make random communities the reference level
allwy$drought <- relevel(allwy$drought, ref = "cntl") #make random communities the reference level


# hist(allwy$dist)
# hist(sqrt(allwy$dist))
#allwy$dist_tran <- sqrt(allwy$dist)

#w/out transfom
summary(t <- aov(dist~trt*drought, allwy))
tuktest <- TukeyHSD(t)
multcompView::multcompLetters4(t,tuktest)
ggplot(allwy, aes(y=dist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Distance from target (normalized)")+ #, fill="drought treatment")+
  theme_classic()

letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought:year'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\2", rownames(letterstest)))
letterstest$year <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\3", rownames(letterstest)))
test <- allwy %>% group_by(year, drought, trt) %>% summarise(yposition = quantile(dist,.8))
test <- merge(letterstest,test, by = c("year", "drought", "trt"))
test2 <- merge(test,allwy, by = c("year", "drought", "trt"), all=T)


ggplot(test2, aes(y=dist, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Distance from target (normalized)")+ #, fill="drought treatment")+
  theme_classic()


# testing 2/6
letterstest <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt:drought'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+)$", "\\2", rownames(letterstest)))
test <- allwy %>% group_by(drought, trt) %>% summarise(yposition = quantile(dist,.8))
test <- merge(letterstest,test, by = c("drought", "trt"))
test2 <- merge(test,allwy, by = c("drought", "trt"), all=T)

droughtcolswy <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color


#export for short report
tiff("figures/cwm wy/distanceplot_wy.tiff", res=400, height = 4,width =8.5, "in",compression = "lzw")
ggplot(test2, aes(y=dist, x=trt, fill=drought))+
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
dev.off()

#combined with distance by year on other r script
alldistwy <- ggplot(test2, aes(y=dist, x=trt, fill=drought))+
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
