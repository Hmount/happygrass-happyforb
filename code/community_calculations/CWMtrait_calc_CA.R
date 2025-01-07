#### Calculate the CWM traits for each plot in each year based on the 
#### taxonomic composition for California site. 
#### Note: in 2022 & 2023 plots are removed that have >20% cover of non-native
#### volunteer species for which we do not have traits data. 

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

#### CA ####
## load in composition data, clean and modify columns as usual
comp.ca <- read.csv("data/raw_cover/Species_Composition_allyears.csv") #read in California comp data
comp.ca <- comp.ca[,c(1:57)]
#comp.ca <- comp.ca %>% filter(Year == "2022") #keep only 2022 data right now
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)


# CSV of species-trait combinations (for OG 25)
traits.ca <- read.csv("data/annualgrass.csv", header=TRUE, row.names=1)
#traits.ca <- traits.ca[traits.ca$use==1,] # subset use=1
#traits$PLSg.m2.mono <- traits$PLSlbs.acre.mono * (453.59237 / 4046.8564224) #convert lb/acre to g/m2
#traits.ca <- traits.ca %>% select(-c(Asat,WUE,RLD,RTD))
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

#filter for species used
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

### pca
pca <- princomp(traits.ca[,2:8], cor=TRUE)
summary(pca)
biplot(pca)
traits.ca$pc1 <- pca$scores[,1]
traits.ca$pc2 <- pca$scores[,2]

#set color scheme
colors = brewer.pal(12,"Paired")
traits.ca$cols = c(colorRampPalette(colors)(nrow(traits.ca)))

### Subset monocots and dicots
grams <- subset(traits.ca, graminoid==1)
forbs <- subset(traits.ca, graminoid==0)

# Arrange trait matrix alphabetically
#trait.matrix.ca <- traits.ca[order(traits.ca$Code),]
trait.matrix.ca <- traits.ca[order(rownames(traits.ca)),]
# # remove species not present in other comp dataframe
# trait.matrix.ca <- trait.matrix.ca %>% 
#   filter(Code != "AVBA") %>% 
#   filter(Code != "BRNI") %>%
#   filter(Code != "CAME") %>%
#   #filter(Code != "FEMY") %>% 
#   filter(Code != "LUAL") %>%
#   filter(Code != "MASA") %>% 
#   filter(Code != "PEHE") %>% 
#   filter(Code != "SACO") 

#Select columns of interest, make into matrix
test <- data.frame(trait.matrix.ca[,2:8], row.names = trait.matrix.ca[,1])
rownames(test)
trait.matrix.ca <- as.matrix(test)
trait.matrix.ca <- trait.matrix.ca[order(rownames(trait.matrix.ca)),]

## pretreatment/ seeding probability 
# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
preds.ca <- read.csv("data/calgrass.allplot.assemblages.csv") #data
preds.ca$trt.b <- paste(preds.ca$trt, preds.ca$block)
preds.ca <- preds.ca %>% arrange(trt.b,block)
comms_p.ca <- labdsv::matrify(data.frame(preds.ca$trt.b,preds.ca$species,preds.ca$prob))
comms_p.ca <- comms_p.ca[,order(colnames(comms_p.ca))]
#comms_p.ca <- comms_p.ca %>% select(-CAME)
trait.matrix.ca.pred <- traits.ca[order(rownames(traits.ca)),]
# remove species not present in other comp dataframe
# trait.matrix.ca.pred <- trait.matrix.ca.pred %>% 
#   filter(Code != "AVBA") %>% 
#   filter(Code != "BRNI") %>%
#   filter(Code != "LUAL") %>%
#   filter(Code != "MASA") %>% 
#   filter(Code != "PEHE") %>% 
#   filter(Code != "SACO")
test <- data.frame(trait.matrix.ca.pred[,c(2:4,7,9,10:12)], row.names = trait.matrix.ca.pred[,1])
rownames(test)
trait.matrix.ca.pred <- as.matrix(test)
trait.matrix.ca.pred <- trait.matrix.ca.pred[order(rownames(trait.matrix.ca.pred)),]
# remove22.pred <- c("FEMY","FEPE","BRMA")
# trait.matrix.ca.pred <- trait.matrix.ca.pred[!row.names(trait.matrix.ca.pred)%in%remove22.pred,]
cwm_p.ca <- FD::functcomp(as.matrix(trait.matrix.ca.pred), as.matrix(comms_p.ca), bin.num=c("graminoid"))
#Define block by extracting the numeric from the cwm rownames
cwm_p.ca$block <- as.factor(sub("^[a-z]+ (\\d+)$", "\\1", rownames(cwm_p.ca)))
#Define seeding trt by extracting the letters from the cwm rownames
cwm_p.ca$trt <- as.factor(sub("(^[a-z]+) \\d+$", "\\1", rownames(cwm_p.ca)))
cwm_p.ca <- cwm_p.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand 
# Define drought treatment at block level
block.water <- comp.ca %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
block.water <- block.water[c(1:210),] #repeating 3 times(?) 
cwm_p.ca <- merge(cwm_p.ca, block.water, by = c("block","trt"), all.x=T) #merge


## 2021
#Arrange cover estimates for field data by year
hpf21 <- comp.ca %>% filter(Year == "2021") %>% arrange(trt,block) 
hpf21 <- hpf21 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

#select cover columns only for matrix
comms21.ca <- hpf21[,c(4,19:57)]

# make comms21.ca the long to add 4-letter codes
comms21.ca.long <- comms21.ca %>% pivot_longer(cols = c(2:40),names_to = "X6letter", values_to = "cover")
comms21.ca.long <- comms21.ca.long %>% mutate(sppcodes = paste0(substr(X6letter, 1, 2), substr(X6letter, 4, 5))) # make 4 letter code to get cwm's

# return to wide matrix
# Reshape field data into a matrix and arrange alphabetically
comms21.ca <- labdsv::matrify(data.frame(comms21.ca.long$trt.b,comms21.ca.long$sppcodes,comms21.ca.long$cover))
comms21.ca <- comms21.ca[,order(colnames(comms21.ca))]
colnames(comms21.ca)[colnames(comms21.ca) == "CACI"] <- "CAME" #change name of CACI to CAME

# remove columns from comms (communties) that are not present in trait.matrix (can we get these? too small proportion?)
comms21.ca <- comms21.ca %>% select(-c("BRHO", "BRDI", "HOMU", "BRMA", "FEPE"))
#remove species from matrix with NA in all communities this year (not present yet) 
#trait.matrix.ca21 <- trait.matrix.ca[!row.names(trait.matrix.ca)=="FEMY",]
# rownames(trait.matrix.ca22)
# colnames(comms21.ca)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm21.ca <- FD::functcomp(as.matrix(trait.matrix.ca), as.matrix(comms21.ca), bin.num=c("graminoid"))

##plots with all 0 need to be removed to run code, but this will change when we get weed trait data
cwm21.ca <- cwm21.ca[rowSums(is.na(cwm21.ca)) != ncol(cwm21.ca), ]
#right now, changing cwm data should be fine, but if these data are needed later, use removed22.wy to progogate figures. 
#removed21.ca <- cwm21.ca[rowSums(is.na(cwm21.ca)) == ncol(cwm21.ca), ] # 2fd, 2ir, 2r, 1dt (see removed plots)
# Building out treatment identification which was absent from our original csv.
N <- 53
#groups.ca <- c(rep("dt",N), rep("fd",N-1), rep("ir",N-1), rep("rand",N))
groups.ca <-c(rep("dt",N-2), rep("fd",N-3), rep("ir",N-3), rep("rand",N-2)) # groups removing uncalculatable rows

#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms21 <- vegan::metaMDS(cwm21.ca, distance="euclidean")
nms21
plot(nms21, type="t", main="euclidean")
ordiellipse(nms21, groups.ca, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms21,cwm21.ca)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm21.ca$trt <- factor(groups.ca)
# cwm_p$trt <- factor(groups)
cwm21.ca$trt <- factor(cwm21.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm21.ca$Treatments <- cwm21.ca$trt
levels(cwm21.ca$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#remove species not occurring in any community this year
sort(colSums(comms21.ca)) #find species
comms21.ca <- comms21.ca %>% select(-c(ELCO,LAPL,MEIM))
remove21 <- c("ELCO", "LAPL", "MEIM")
trait.matrix.ca21 <- trait.matrix.ca[!row.names(trait.matrix.ca)%in%remove21,]

#remove communities that do not have any species data in 2021
comms21.ca <- comms21.ca %>%
  filter(rowSums(comms21.ca) != 0)

#Calculating functional diversity of each trait and extracting RaoQ on L115
raoq <- FD::dbFD(as.matrix(trait.matrix.ca21), as.matrix(comms21.ca))
leafn <- FD::dbFD(as.matrix(trait.matrix.ca21[,"N"]), as.matrix(comms21.ca))
lma <- FD::dbFD(as.matrix(trait.matrix.ca21[,"LMA"]), as.matrix(comms21.ca))
seedmass <- FD::dbFD(as.matrix(trait.matrix.ca21[,"seed.mass"]), as.matrix(comms21.ca))
srl <- FD::dbFD(as.matrix(trait.matrix.ca21[,"SRL"]), as.matrix(comms21.ca))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.ca21[,"Rdiam"]), as.matrix(comms21.ca))
rmf <- FD::dbFD(as.matrix(trait.matrix.ca21[,"RMF"]), as.matrix(comms21.ca))
cwm_roaq21.ca <- data.frame(leafn=leafn$RaoQ,
                       lma=lma$RaoQ,
                       seedmass=seedmass$RaoQ,
                       srl=srl$RaoQ,
                       rootdiam=rootdiam$RaoQ,
                       rmf=rmf$RaoQ,
                       full=raoq$RaoQ,
                       trt=factor(groups.ca))
cwm_roaq21.ca$trt <- factor(cwm_roaq21.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm21.ca$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm21.ca)))
# Define drought treatment at block level
block.water <- hpf21 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cwm21.ca <- merge(cwm21.ca, block.water, by = c("block","trt"), all.x=T) #merge



## 2022
#Arrange cover estimates for field data by year
hpf22 <- comp.ca %>% filter(Year == "2022") %>% arrange(trt,block) 
hpf22 <- hpf22 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

#select cover columns only for matrix
comms22.ca <- hpf22[,c(4,19:57)]

# make comms22.ca the long to add 4-letter codes
comms22.ca.long <- comms22.ca %>% pivot_longer(cols = c(2:40),names_to = "X6letter", values_to = "cover")
comms22.ca.long <- comms22.ca.long %>% mutate(sppcodes = paste0(substr(X6letter, 1, 2), substr(X6letter, 4, 5))) # make 4 letter code to get cwm's

# return to wide matrix
# Reshape field data into a matrix and arrange alphabetically
comms22.ca <- labdsv::matrify(data.frame(comms22.ca.long$trt.b,comms22.ca.long$sppcodes,comms22.ca.long$cover))
comms22.ca <- comms22.ca[,order(colnames(comms22.ca))]
colnames(comms22.ca)[colnames(comms22.ca) == "CACI"] <- "CAME" #change name of CACI to CAME

# remove columns from comms (communties) that are not present in trait.matrix (can we get these? too small proportion?)
comms22.ca <- comms22.ca %>% select(-c("BRMA", "BRDI","FEPE"))
#remove species from matrix with NA in all communities this year (not present yet) 
#trait.matrix.ca22 <- trait.matrix.ca[!row.names(trait.matrix.ca)=="FEMY",]
 rownames(trait.matrix.ca)
 colnames(comms22.ca)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm22.ca <- FD::functcomp(as.matrix(trait.matrix.ca), as.matrix(comms22.ca), bin.num=c("graminoid"))

# Building out treatment identification which was absent from our original csv.
N <- 53
# plots with all 0 (no spp in community) need to be removed to run nms2,
cwm22.ca <- cwm22.ca[rowSums(is.na(cwm22.ca)) != ncol(cwm22.ca), ] 
#right now, changing cwm22.ca data should be fine, but if these data are needed later, use removed23.wy to propagate figures. 
removed22.ca <- cwm22.ca[rowSums(is.na(cwm22.ca)) == ncol(cwm22.ca), ] # 2fd,3ir,1dt (see removed plots)
groups.ca <-c(rep("dt",N-1), rep("fd",N-3), rep("ir",N-4), rep("rand",N)) # groups removing uncalculatable rows

#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms22.ca <- vegan::metaMDS(cwm22.ca, distance="euclidean")
nms22.ca
plot(nms22.ca, type="t", main="euclidean")
ordiellipse(nms22.ca, groups.ca, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms22.ca,cwm22.ca)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm22.ca$trt <- factor(groups.ca)
# cwm_p$trt <- factor(groups)
cwm22.ca$trt <- factor(cwm22.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm22.ca$Treatments <- cwm22.ca$trt
levels(cwm22.ca$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#remove communities that do not have any species data in 2022
comms22.ca <- comms22.ca %>%
  filter(rowSums(comms22.ca) != 0) #remove comms that sum to 0

#remove species not occurring in any community this year
sort(colSums(comms22.ca)) #find species, 5 spp = 0
comms22.ca <- comms22.ca %>% select(-c(ARPU,KOMA,MURI,TRWI,MEIM))
remove22.ca <- c("ARPU", "KOMA", "MURI", "TRWI", "MEIM")
trait.matrix.ca22 <- trait.matrix.ca[!row.names(trait.matrix.ca)%in%remove22.ca,]

#Calculating functional diversity of each trait and extracting RaoQ on L115
raoq <- FD::dbFD(as.matrix(trait.matrix.ca22), as.matrix(comms22.ca))
leafn <- FD::dbFD(as.matrix(trait.matrix.ca22[,"N"]), as.matrix(comms22.ca))
lma <- FD::dbFD(as.matrix(trait.matrix.ca22[,"LMA"]), as.matrix(comms22.ca))
seedmass <- FD::dbFD(as.matrix(trait.matrix.ca22[,"seed.mass"]), as.matrix(comms22.ca))
srl <- FD::dbFD(as.matrix(trait.matrix.ca22[,"SRL"]), as.matrix(comms22.ca))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.ca22[,"Rdiam"]), as.matrix(comms22.ca))
rmf <- FD::dbFD(as.matrix(trait.matrix.ca22[,"RMF"]), as.matrix(comms22.ca))
cwm_roaq22.ca <- data.frame(leafn=leafn$RaoQ,
                       lma=lma$RaoQ,
                       seedmass=seedmass$RaoQ,
                       srl=srl$RaoQ,
                       rootdiam=rootdiam$RaoQ,
                       rmf=rmf$RaoQ,
                       full=raoq$RaoQ,
                       trt=factor(groups.ca))
cwm_roaq22.ca$trt <- factor(cwm_roaq22.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm22.ca$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm22.ca)))
# Define drought treatment at block level
block.water <- hpf22 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cwm22.ca <- merge(cwm22.ca, block.water, by = c("block","trt"), all.x=T) #merge



## 2023
#Arrange cover estimates for field data by year
hpf23 <- comp.ca %>% filter(Year == "2023") %>% arrange(trt,block) 
hpf23 <- hpf23 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

#select cover columns only for matrix
comms23.ca <- hpf23[,c(4,19:57)]

# make comms23.ca the long to add 4-letter codes
comms23.ca.long <- comms23.ca %>% pivot_longer(cols = c(2:40),names_to = "X6letter", values_to = "cover")
comms23.ca.long <- comms23.ca.long %>% mutate(sppcodes = paste0(substr(X6letter, 1, 2), substr(X6letter, 4, 5))) # make 4 letter code to get cwm's

# return to wide matrix
# Reshape field data into a matrix and arrange alphabetically
comms23.ca <- labdsv::matrify(data.frame(comms23.ca.long$trt.b,comms23.ca.long$sppcodes,comms23.ca.long$cover))
comms23.ca <- comms23.ca[,order(colnames(comms23.ca))]
colnames(comms23.ca)[colnames(comms23.ca) == "CACI"] <- "CAME" #change name of CACI to CAME

# remove columns not present in communities this year
comms23.ca <- comms23.ca %>% select(-c("BRMA", "BRDI", "BRHO", "HOMU", "POMO","FEMY","FEPE"))
rownames(trait.matrix.ca)
colnames(comms23.ca)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm23.ca <- FD::functcomp(as.matrix(trait.matrix.ca), as.matrix(comms23.ca), bin.num=c("graminoid"))

# Building out treatment identification which was absent from our original csv.
N <- 53

##plots with all 0 (no spp in community) need to be removed to run nms2,
cwm23.ca <- cwm23.ca[rowSums(is.na(cwm23.ca)) != ncol(cwm23.ca), ] 
#right now, changing cwm data should be fine, but if these data are needed later, use removed23.wy to propagate figures. 
#removed23.ca <- cwm23.ca[rowSums(is.na(cwm23.ca)) == ncol(cwm23.ca), ] # fd.36 (see removed plot)
groups.ca <-c(rep("dt",N-1), rep("fd",N-3), rep("ir",N-3), rep("rand",N-1)) # groups removing uncalculatable rows

#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms23.ca <- vegan::metaMDS(cwm23.ca, distance="euclidean")
nms23.ca
plot(nms23.ca, type="t", main="euclidean")
ordiellipse(nms23.ca, groups.ca, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms23.ca,cwm23.ca)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm23.ca$trt <- factor(groups.ca)
# cwm_p$trt <- factor(groups)
cwm23.ca$trt <- factor(cwm23.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm23.ca$Treatments <- cwm23.ca$trt
levels(cwm23.ca$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#remove communities that do not have any species data in 2023
comms23.ca <- comms23.ca %>%
  filter(rowSums(comms23.ca) != 0)

#remove species not occurring in any community this year
sort(colSums(comms23.ca)) #find species, 5 spp = 0
comms23.ca <- comms23.ca %>% select(-c(ARPU,KOMA,LAPL,MURI,TRWI,MEIM))
remove23.ca <- c("ARPU", "KOMA","LAPL", "MURI", "TRWI", "MEIM")
trait.matrix.ca23 <- trait.matrix.ca[!row.names(trait.matrix.ca)%in%remove23.ca,]

#Calculating functional diversity of each trait and extracting RaoQ on L115
raoq <- FD::dbFD(as.matrix(trait.matrix.ca23), as.matrix(comms23.ca))
leafn <- FD::dbFD(as.matrix(trait.matrix.ca23[,"N"]), as.matrix(comms23.ca))
lma <- FD::dbFD(as.matrix(trait.matrix.ca23[,"LMA"]), as.matrix(comms23.ca))
seedmass <- FD::dbFD(as.matrix(trait.matrix.ca23[,"seed.mass"]), as.matrix(comms23.ca))
srl <- FD::dbFD(as.matrix(trait.matrix.ca23[,"SRL"]), as.matrix(comms23.ca))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.ca23[,"Rdiam"]), as.matrix(comms23.ca))
rmf <- FD::dbFD(as.matrix(trait.matrix.ca23[,"RMF"]), as.matrix(comms23.ca))
cwm_roaq23.ca <- data.frame(leafn=leafn$RaoQ,
                       lma=lma$RaoQ,
                       seedmass=seedmass$RaoQ,
                       srl=srl$RaoQ,
                       rootdiam=rootdiam$RaoQ,
                       rmf=rmf$RaoQ,
                       full=raoq$RaoQ,
                       trt=factor(groups.ca))
cwm_roaq23.ca$trt <- factor(cwm_roaq23.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm23.ca$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm23.ca)))
# Define drought treatment at block level
block.water <- hpf23 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cwm23.ca <- merge(cwm23.ca, block.water, by = c("block","trt"), all.x=T) #merge


## merge cwm's together for storing
cpre <- cwm_p.ca %>% mutate(year = "0") 
cpre$trt <- as.factor(as.character(cpre$trt))
rownames(cpre) <- NULL
c21 <- cwm21.ca %>% mutate(year = "2021") #add year
c21$trt <- as.factor(as.character(c21$trt))
rownames(c21) <- NULL
c22 <- cwm22.ca %>% mutate(year = "2022") #add year
c22$trt <- as.factor(as.character(c22$trt))
rownames(c22) <- NULL
c23 <- cwm23.ca %>% mutate(year = "2023") #add year
c23$trt <- as.factor(as.character(c23$trt))
rownames(c23) <- NULL

cwm.ca <- bind_rows(cpre, c21) #bind 1st
cwm.ca <- bind_rows(cwm.ca, c22) #bind again
cwm.ca <- bind_rows(cwm.ca, c23) #bind again

### how many WY 2022+2023 data need to be dropped from CWM calculations 
### important column unless we get trait data
subca <- comp.ca %>% #[,-c(18:56)]
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>%
  filter(propnative < 80)
table(subca$year)
table(comp.ca$year)
(13+22+62)/(210*3)*100 # only 15% total observation to remove
(62)/(210)*100 # 30% total observation to remove for inv. models
subca <- subca[,c(1,4,5,59)]
subca$trt <- tolower(subca$trt) #make lowercase to match
subca$block <- as.factor(subca$block) #make lowercase to match
subca$Year <- as.factor(subca$Year) #make lowercase to match

#merge
test <- merge(cwm.ca,subca, by.x=c("block","trt","year"), by.y=c("block","trt","Year"), all.x=T)

#save 
write.csv(cwm.ca, "data/cwm_ca.csv", row.names = F)

## merge cwm's together for storing
cfd1 <- cwm_roaq21.ca %>% mutate(year = "2021") #add year
cfd1$block <- as.numeric(sub(".*\\.(\\d+)$", "\\1", rownames(cfd1))) #Define block by extracting the numeric from the cwm rownames
# Define drought treatment at block level
block.water <- hpf22 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cfd1 <- merge(cfd1, block.water, by = c("block","trt"), all.x=T) #merge

cfd2 <- cwm_roaq22.ca %>% mutate(year = "2022") #add year
cfd2$block <- as.numeric(sub(".*\\.(\\d+)$", "\\1", rownames(cfd2))) #Define block by extracting the numeric from the cwm rownames
cfd2 <- merge(cfd2, block.water, by = c("block","trt"), all.x=T) #merge

cfd3 <- cwm_roaq23.ca %>% mutate(year = "2023") #add year
cfd3$block <- as.numeric(sub(".*\\.(\\d+)$", "\\1", rownames(cfd3))) #Define block by extracting the numeric from the cwm rownames
cfd3 <- merge(cfd3, block.water, by = c("block","trt"), all.x=T) #merge

cwm.cafd <- bind_rows(cfd1, cfd2) #bind 1st
cwm.cafd <- bind_rows(cwm.cafd, cfd3) #bind again
#save 
write.csv(cwm.cafd, "data/cwm_raoq_ca.csv", row.names = F)
