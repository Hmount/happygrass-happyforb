#### CWM calculation, Euclidean distance calculation, and distance figure for 2020 
#### (pre-treatment) in WY.
#### We were unable to collect this data in California due to the pandemic
#### response to reviewer comments, but not included in manuscipt or appendix

## packages
library(tidyverse)
library(FD)
library(mice)
library(RColorBrewer)
library(rlist)
library(labdsv)
library(vegan)
library(sads)
library(erer)

### CWM first:
## load in composition data, clean and modify columns as usual
comp.wy <- read.csv("data/comp_wy_plot.csv") #Wyoming species comp data PLOT level
comp.wy <- comp.wy %>% filter(year == "2020") #remove 2020 data 
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy <- comp.wy %>% unite(plot, c(block, trt), sep = ".", remove=F) # make unique plot variable

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

#set color scheme
colors = brewer.pal(12,"Paired")
traits.wy$cols = c(colorRampPalette(colors)(nrow(traits.wy)))

#this code is from calculating year 0 in CWMtrait_calc_WY(plot).R script
#og25 <- unique(colnames(comms_p.wy)) #make list of 25 native species

### Subset monocots and dicots
grams <- subset(traits.wy, graminoid==1)
forbs <- subset(traits.wy, graminoid==0)

# Arrange trait matrix alphabetically
trait.matrix.wy <- traits.wy[order(rownames(traits.wy)),]
#Select colons of interest
trait.matrix.wy <- as.matrix(trait.matrix.wy[,2:11])

#Arrange cover estimates for field data by year and remove weed species for now
hpg20 <- comp.wy %>% arrange(trt,block)
hpg20 <- hpg20 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
needcols <- c(og25,"trt.b")
comms20.wy <- hpg20 %>% column_to_rownames("trt.b")
comms20.wy <- comms20.wy %>% select(all_of(og25)) #ONLY NATIVE 25 and grouping columns
comms20.wy <- comms20.wy[,order(colnames(comms20.wy))]
comms20.wy <- replace(comms20.wy, is.na(comms20.wy), 0)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm20.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms20.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64*1
groups.wy <- c(rep("dt",N), rep("fd",N), rep("ir",N), rep("rand",N))

#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms20.wy <- vegan::metaMDS(cwm20.wy, distance="euclidean")
nms20.wy
plot(nms20.wy, type="t", main="euclidean")
ordiellipse(nms20.wy, groups.wy, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms20.wy,cwm20.wy)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm20.wy$trt <- factor(groups.wy)
# cwm_p$trt <- factor(groups)
cwm20.wy$trt <- factor(cwm20.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm20.wy$Treatments <- cwm20.wy$trt
levels(cwm20.wy$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#remove ARPU since it was never present and dbFD won't run
comms20.wy <- comms20.wy[,colnames(comms20.wy) != "ARPU"]
trait.matrix.wy <- trait.matrix.wy[rownames(trait.matrix.wy) != "ARPU", ]
#Calculating functional diversity of each trait and extracting RaoQ (RoaQ is not scaled here, needs to be done when using these columns)
raoq <- FD::dbFD(as.matrix(trait.matrix.wy), as.matrix(comms20.wy))
leafn <- FD::dbFD(as.matrix(trait.matrix.wy[,"leafn"]), as.matrix(comms20.wy))
lop <- FD::dbFD(as.matrix(trait.matrix.wy[,"lop"]), as.matrix(comms20.wy))
ldmc <- FD::dbFD(as.matrix(trait.matrix.wy[,"ldmc"]), as.matrix(comms20.wy))
srl <- FD::dbFD(as.matrix(trait.matrix.wy[,"srl"]), as.matrix(comms20.wy))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.wy[,"rootdiam"]), as.matrix(comms20.wy))
veg <- FD::dbFD(as.matrix(trait.matrix.wy[,"veg"]), as.matrix(comms20.wy))
cwm_roaq20.wy <- data.frame(leafn=leafn$RaoQ,
                            lop=lop$RaoQ,
                            ldmc=ldmc$RaoQ,
                            srl=srl$RaoQ,
                            rootdiam=rootdiam$RaoQ,
                            veg=veg$RaoQ,
                            full=raoq$RaoQ,
                            trt=factor(groups.wy))
cwm_roaq20.wy$trt <- factor(cwm_roaq20.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm20.wy$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm20.wy)))
cwm_roaq20.wy$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm_roaq20.wy)))
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))


## IR
distir.wy <- cwm20.wy %>% select(c(block,trt,leafn,srl))
distir.wy <- merge(distir.wy,cwm_roaq20.wy[,c(6,8:9)], by=c("trt", "block"),all.x = T)
distir.wy$veg <- normalize(distir.wy$veg) #function in distance_models_figures.R script
# within each year subset the data
irdist20 <- distir.wy %>% unite(trt.b.y, c(trt, block), sep = ".", remove=T)
# attach row of targets
irdist20 <- irdist20 %>% add_row(trt.b.y = "target",
                                 leafn = quantile(irdist20$leafn,.05),
                                 srl = quantile(irdist20$srl,.95),
                                 veg = quantile(irdist20$veg,.95))
irdist20 <- irdist20 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
irdistmat20 <- vegdist(as.matrix(irdist20),method = "euclidean", upper=T)#,diag=T)
irdistmat20 <-as.matrix(irdistmat20)
# save only pairwise between target
irdistances20 <- as.data.frame(irdistmat20["target",])
colnames(irdistances20) <- "distir"

## DT
distdt.wy <- cwm20.wy %>% select(c(block,trt,ldmc,lop))
distdt.wy <- merge(distdt.wy,cwm_roaq20.wy[,c(5,8:9)], all.x = T)
# within each year subset the data
dtdist20 <- distdt.wy %>% unite(trt.b.y, c(trt, block), sep = ".", remove=T)
# attach row of targets
dtdist20 <- dtdist20 %>% add_row(trt.b.y = "target",
                                 ldmc = quantile(dtdist20$ldmc,.95),
                                 lop = quantile(dtdist20$lop,.05),
                                 rootdiam = quantile(dtdist20$rootdiam,.95))
dtdist20 <- dtdist20 %>% column_to_rownames("trt.b.y")
#run dist or vegdist
dtdistmat20 <- vegdist(as.matrix(dtdist20),method = "euclidean", upper=T)#,diag=T)
dtdistmat20 <-as.matrix(dtdistmat20)
# save only padtwise between target
dtdistances20 <- as.data.frame(dtdistmat20["target",])
colnames(dtdistances20) <- "distdt"

#FD
fddist <- cwm_roaq20.wy[,c(7,8:9)]# %>% filter(trt=="fd")
#FD target (shifting annually)
quantile(normalize(fddist$full),.99)
#using max(), gives 1 in all years due to normalization, 99th quantile shows variation
subdatFD <- fddist %>% summarize(target = quantile(normalize(full),.99))
#quantile is higher when considering less plots (i.e. only FD) and more representative of the
# comparison between FD and target FD so using subsetted data to calculate dist
# fddist <- distfd.wy %>%
#   unite(trt.b.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
fddist$full <- normalize(fddist$full)
##2020
# within each year subset the data
fddist20 <- fddist %>% unite(trt.b.y, c(trt, block), sep = ".", remove=T)
# attach row of targets
fddist20 <- fddist20 %>% add_row(trt.b.y = "target",
                                 full = quantile(fddist20$full,.99)) #this should be by year tho
fddist20 <- fddist20 %>% select(-trt.b.y)
#run dist or vegdist
fddistmat <- vegdist(as.matrix(fddist20),method = "euclidean", upper=T)#,diag=T)
fddistmat <-as.matrix(fddistmat)
# save only pairwise between target
fddistances20 <- as.data.frame(fddistmat["...257",]) #target
colnames(fddistances20) <- "distfd"

rdist <- cwm20.wy %>% arrange(block)
rdist$trt <- factor(rdist$trt, ordered = FALSE)
#rdist0 <- rdist %>% filter(year=="0") %>% arrange(as.numeric(block)) #from distance_calculatino_WY.R

#2020
rdist20 <- bind_rows(rdist0,rdist)
rdist20 <- rdist20 %>% select(-c(drought, Treatments))
rdist20 <- rdist20 %>% unite(trt.b.y, c(trt,block,year), sep = ".", remove=T) # make unique plot variable
rdist20 <- rdist20 %>% select(-trt.b.y)
randdistmat <- vegdist(as.matrix(rdist20),method = "euclidean", upper = T)
rdist20 <-as.matrix(randdistmat)[,-c(1:256)]
rdist20 <-rdist20[c(1:256),]
rdist20<- data.frame(
  dist=diag(as.matrix(rdist20)),
  id=colnames(rdist20))
rdistances <- rdist20
colnames(rdistances) <- c("distr","trt.b.y")
rdistances$trt.b.y <- as.factor(gsub("^r", "rand", rdistances$trt.b.y))

irdistances20 <- irdistances20 %>% rownames_to_column("trt.b.y")
dtdistances20 <- dtdistances20 %>% rownames_to_column("trt.b.y")
fddistances20 <- fddistances20 %>% rownames_to_column("trt.b.y")
fddistances20$trt.b.y <- as.factor(gsub("^r", "rand", fddistances20$trt.b.y))

#### combine using min/max in DT and IR
wydist20 <- merge(dtdistances20,irdistances20)
wydist20 <- merge(wydist20,fddistances20)
wydist20 <- merge(wydist20,rdistances)


### model/ figure:
## seperate plot ID column into trt, block, and year
wydist <- wydist20 %>% 
  separate(trt.b.y, into = c("trt", "block"), sep = "\\.")
wydist$block <- as.numeric(wydist$block)
## set reference levels for modelling
wydist$trt <- as.factor(wydist$trt) #must be factor
#wydist$trt <- relevel(wydist$trt, ref = "rand") #make random communities the reference level

## Drought Tolerant
## are drought tolerant plots significantly closer to our DT target than random?
wydist.dt <- wydist %>% filter(trt=="dt"|trt=="rand")# subset for only DT and RC communities
summary(t <- aov(distdt~trt, wydist.dt)) #run model 
#create letters for plotting:
tuktest <- TukeyHSD(t)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
dttemp <- wydist.dt %>% group_by(trt) %>% summarise(yposition = quantile(distdt,.8))
dttemp <- merge(letters,dttemp, by =  "trt")
dttemp <- merge(dttemp,wydist.dt, by = "trt", all=T)
#plot:
distdtwy <- ggplot(dttemp, aes(y=distdt, x=trt))+
  geom_boxplot()+
  scale_x_discrete(labels = c("DT", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            hjust = 1.5,
            size=3) +
  labs(x=" ",y="... DT target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

## Invasion resistant 
## are invasion resistant plots significantly closer to our IR target than random?
wydist.ir <- wydist %>% filter(trt=="ir"|trt=="rand")# subset for only IR and RC communities
summary(t <- aov(distir~trt, wydist.ir)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t)
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
irtemp <- wydist.ir %>% group_by(trt) %>% summarise(yposition = quantile(distir,.8))
irtemp <- merge(letters,irtemp, by = "trt")
irtemp <- merge(irtemp,wydist.ir, by = "trt", all=T)
#plot:
distirwy <- ggplot(irtemp, aes(y=distir, x=trt))+
  geom_boxplot()+
  scale_x_discrete(labels = c("IR", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.2), 
            vjust = -0.5,
            hjust = 1.5,
            size=3) +
  labs(x=" ",y="... IR target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

## Functionally Diverse
## are these plots significantly more functionally diverse (Rao) than random?
wydist.fd <- wydist %>% filter(trt=="fd"|trt=="rand")# subset for only FD and RC communities
summary(t <- aov(distfd~trt, wydist.fd)) #run model *reported in ms*
#create letters for plotting:
tuktest <- TukeyHSD(t) #run post-hoc for letters
letters <- data.frame(multcompView::multcompLetters4(t,tuktest)$'trt'['Letters'])
letters$trt <- as.factor(sub("([a-z]+):([a-z]+)$", "\\1", rownames(letters)))
fdtemp <- wydist.fd %>% group_by(trt) %>% summarise(yposition = quantile(distfd,.8))
fdtemp <- merge(letters,fdtemp, by = "trt")
fdtemp <- merge(fdtemp,wydist.fd, by = "trt", all=T)
#plot:
distfdwy <- ggplot(fdtemp, aes(y=distfd, x=trt))+
  geom_boxplot()+
  scale_x_discrete(labels = c("FD", "RC"))+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 1.45), 
            vjust = -0.5,
            hjust = 1.5,
            size=3) +
  labs(x=" ",y="... FD target")+ #, fill="drought treatment")+
  theme_ggeffects()+
  theme(legend.position = "none")+
  ylim(0,4)

wydistplots <- ggarrange(distdtwy, distirwy, distfdwy, ncol=3, nrow=1, 
                         common.legend = T, legend = "none", 
                         labels = c("a","b","c"), hjust=c(-4,-4,-4))
wydistplots2 <- annotate_figure(wydistplots, left=text_grob("Euclidean distance to...", rot=90))
