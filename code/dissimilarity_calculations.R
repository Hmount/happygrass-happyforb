# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

library(tidyverse)
dat <- read.csv("data/cwm_wy.csv")
dat$year <- as.factor(dat$year)
dat$block <- as.factor(dat$block)
dat$trt <- as.factor(dat$trt)
dat$subplot <- as.factor(dat$subplot)
FDdat <- read.csv("data/cwm_raoq_wy.csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)
FDdat$subplot <- as.factor(FDdat$subplot)
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

## Calculating dissimilarity/ distance (euclidian) from targets
## select relevent CWM traits per DT or IR, rows are communities
pred <- dat %>% filter(year=="0")#filter preds
dat <- dat %>% filter(year!="0")#filter preds
FDpred <- FDdat %>% filter(year=="0")#filter preds
FDdat <- FDdat %>% filter(year!="0")#filter preds


#subset for just seeding trt of interest
## IR
distir.wy <- dat %>% select(c(block,trt,subplot,year,leafn,srl))
distir.wy <- merge(distir.wy,FDdat[,c(2,4:6,12)], all.x = T)
# define trait targets
quantile(traits.wy$leafn,.25) #IR
# -0.565257 
quantile(traits.wy$srl,.7557) #IR
# 0.6536232 
irdist <- distir.wy %>% filter(trt=="ir")
irdist <- irdist %>% 
  unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
irdist$veg <- normalize(irdist$veg)
# make row of targets
irdist <- irdist %>% add_row(trt.b.sub.y = "target", 
                           leafn = quantile(traits.wy$leafn,.25),
                           srl = quantile(traits.wy$srl,.7557),
                           veg = quantile(normalize(FDdat$veg),.99)) #this should be by year tho
irdist <- irdist %>% column_to_rownames("trt.b.sub.y")
#mnake into matrix
#run dist or vegdist
library(vegan)
irdistmat <- vegdist(as.matrix(irdist),method = "euclidean", upper=T)#,diag=T)
irdistmat <-as.matrix(irdistmat)
# save only pairwise between target 
irdistances <- as.data.frame(irdistmat["target",])
colnames(irdistances) <- "ir_dist"

## DT
distdt.wy <- dat %>% select(c(block,trt,subplot,year,ldmc,lop))
distdt.wy <- merge(distdt.wy,FDdat[,c(1,4:6,12)], all.x = T)
# define trait targets
quantile(traits.wy$ldmc,.75) # DT
# 0.6766338 
quantile(traits.wy$lop,.25) #DT
# -0.7420366
dtdist <- distdt.wy %>% filter(trt=="dt")
dtdist <- dtdist %>% 
  unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
dtdist$rootdiam <- normalize(dtdist$rootdiam)
# make row of targets
dtdist <- dtdist %>% add_row(trt.b.sub.y = "target", 
                             ldmc = quantile(traits.wy$ldmc,.75),
                             lop = quantile(traits.wy$lop,.25),
                             rootdiam = quantile(normalize(FDdat$rootdiam),.99)) #this should be by year tho
dtdist <- dtdist %>% column_to_rownames("trt.b.sub.y")
#mnake into matrix
#run dist or vegdist
dtdistmat <- vegdist(as.matrix(dtdist),method = "euclidean", upper=T)#,diag=T)
dtdistmat <-as.matrix(dtdistmat)
# save only pairwise between target 
dtdistances <- as.data.frame(dtdistmat["target",])
colnames(dtdistances) <- "dt_dist"

## To assess similarity to highest multivariate FD we use upper 95th percentile
## of each years FD as the target to compare to. 
#FD
fddist <- FDdat[,c(3,4:6,12)] %>% filter(trt=="fd")
#FD target (shifting annually)
quantile(normalize(distfd.wy$full),.99) 
#using max(), gives 1 in all years due to normalization, 99th quantile shows variation
# subdatFD <- FDdat %>% group_by(year) %>% summarize(target = quantile(normalize(full),.99))
subdatFD <- fddist %>% group_by(year) %>% summarize(target = quantile(normalize(full),.99))
#quantile is higher when considering less plots (i.e. only FD) and more representative of the 
# comparison between FD and target FD so using subsetted data to calculate dist
# fddist <- distfd.wy %>% 
#   unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
fddist$full <- normalize(fddist$full)
##2021
# within each year subset the data
fddist21 <- fddist %>% filter(year=="2021")
fddist21 <- fddist21 %>% unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T)
# attach row of targets
fddist21 <- fddist21 %>% add_row(trt.b.sub.y = "target",
                                 full = quantile(normalize(fddist21$full),.99)) #this should be by year tho
fddist21 <- fddist21 %>% column_to_rownames("trt.b.sub.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist21),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances21 <- as.data.frame(fddistmat["target",])
colnames(fddistances21) <- "fd_dist"

##2022
# within each year subset the data
fddist22 <- fddist %>% filter(year=="2022")
fddist22 <- fddist22 %>% unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T)
# attach row of targets
fddist22 <- fddist22 %>% add_row(trt.b.sub.y = "target",
                                 full = quantile(normalize(fddist22$full),.99)) #this should be by year tho
fddist22 <- fddist22 %>% column_to_rownames("trt.b.sub.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist22),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances22 <- as.data.frame(fddistmat["target",])
colnames(fddistances22) <- "fd_dist"

##2023
# within each year subset the data
fddist23 <- fddist %>% filter(year=="2023")
fddist23 <- fddist23 %>% unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T)
# attach row of targets
fddist23 <- fddist23 %>% add_row(trt.b.sub.y = "target",
                                 full = quantile(normalize(fddist23$full),.99)) #this should be by year tho
fddist23 <- fddist23 %>% column_to_rownames("trt.b.sub.y")
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist23),method = "euclidean", upper=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target 
fddistances23 <- as.data.frame(fddistmat["target",])
colnames(fddistances23) <- "fd_dist"

#bind dataframes
fddistances <- bind_rows(fddistances21,fddistances22)
fddistances <- bind_rows(fddistances,fddistances23)


## To assess similarity to random so it can be included in models, we consider the distances 
## in community abundance/ composition from our expected probabilities from seeding.
## this essentially uses the random community as a measure of how well be were at making any 
## community we sought out to make by seeding and allows us to continue using random as a control. 
# Rand
compdat <- read.csv("data/comp_wy.csv")
compdat <- compdat[,c(1:4,11:66)]
rdist <- compdat %>% filter(trt=="r") %>%
  filter(year!="2020")
rdist <- rdist %>% #select(-c(Treatments,subplot))
  unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T) # make unique plot variable

# combine the above subsetted data with the targets (partially done below), use jesse code 
# to run vegdist and extract pairwise differences from mstrix. 
### HERE ###

# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
preds.wy <- read.csv("data/allplot.assemblages.csv") #data
preds.wy <- preds.wy %>% arrange(trt,block)
comms_p.wy <- labdsv::matrify(data.frame(preds.wy$trt,preds.wy$species,preds.wy$prob*100))
comms_p.wy <- comms_p.wy[,order(colnames(comms_p.wy))]
# comms_p.wy$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(comms_p.wy)))
# #Define trt by extracting the subplot from the cwm rownames
# comms_p.wy$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(comms_p.wy))))
# Define drought treatment at block level
# covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
# comms_p.wy <- comms_p.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
#                                                     !block %in% covered ~ "cntl")) 


#2021
# within each year subset the data
fddist21 <- fddist %>% filter(year=="2021")
fddist21 <- fddist21 %>% unite(trt.b.sub.y, c(trt, block, year), sep = ".", remove=T)

randdistmat <- vegdist(as.matrix(rdist),method = "euclidean", upper=T)

ran <- data.frame(rbind(
  #   cbind(str_c(hpg$trt,"E"),hpg$species,hpg$prob*100), # expected coverage from seeded probability
  #   cbind(str_c(hpg$trt,"A"),hpg$species,hpg$cov)))     # actual coverage (was absolute, I have means)
  # ran[,3] <- as.numeric(ran[,3])
  # comms_ran <- labdsv::matrify(ran)
  # comms_ran <- comms_ran[,order(colnames(comms_ran))]
  # comms_ran <- comms_ran %>% select(-BG)
  # 
  # ran_exp <- vegdist(comms_ran[385:512,],upper = T)
  # # Extract all values 1 position below diag and extract every second value
  # ran_dist <- data.frame(
  #   dist=diag(as.matrix(ran_exp)[, -1])[c(TRUE, FALSE)],
  #   id=rownames(comms_ran[385:512,])[c(TRUE, FALSE)])
  # ran_dist$id <- str_remove(ran_dist$id,"A")
  # rownames(ran_dist) <- ran_dist$id
  
  # Random distances are bray curtis distances from the expected targets (prob)
  #the target set out compositionally
  cwm_rand <- left_join(
    cwm %>% filter(trt == "rand") %>% 
      mutate(id = rownames(cwm %>% filter(trt == "rand"))),
    ran_dist, by="id") %>% #make ran_dist not just 2020
    select(-dist_target) %>% 
    rename(dist_target = dist) 
  
  rownames(cwm_rand) <- cwm_rand$id
  cwm_rand <- cwm_rand %>% select(-id)
