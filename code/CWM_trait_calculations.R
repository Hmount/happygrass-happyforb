#### CWM traits calculations, each site and year are seperatley calculated
#### WY first. 2021 ready, but 2022 & 2023 need trait data for volunteer invaders still

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

## functions for traits plots 
CWM_trait_plot <- function(data, trait, points_data) {
  ggplot(data, aes_string("trt", trait)) +
    geom_violin(aes(fill = trt), width = 1, trim = FALSE) +
    geom_boxplot(width = 0.25) +
    geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.3) +
    geom_point(
      aes_string(y = paste0("quantile(", trait, ",.25)"), x = "ir"),
      data = points_data,
      col = "red",
      shape = 18,
      size = 3.5
    ) +
    scale_fill_viridis_d(option = "D", begin = 0.1, end = 1, alpha = 0.7) +
    theme_classic() +
    ylab(trait) +
    xlab("") +
    theme(legend.position = "none")
}

## functions for traits functional diversity (RoaQ).
CWM_trait_FD_plot <- function(data, trait, ylab) {
  ggplot(data, aes_string("trt", trait)) +
    geom_violin(aes(fill = trt), width = 1, trim = FALSE) +
    geom_boxplot(width = 0.25) +
    geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.3) +
    scale_fill_viridis_d(option = "D", begin = 0.1, end = 1, alpha = 0.7) +
    theme_classic() +
    ylab(paste("FD (RaoQ) of ", ylab)) +
    xlab("") +
    ylim(0,max(trait+.5)) +
    theme(legend.position = "none")
}

## functions for trait plots 
CWM_trait_plot <- function(data, trait, ylab) {
  ggplot(data, aes_string("trt", trait)) +
    geom_violin(aes(fill = trt), width = 1, trim = FALSE) +
    geom_boxplot(width = 0.25) +
    geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.3) +
    scale_fill_viridis_d(option = "D", begin = 0.1, end = 1, alpha = 0.7) +
    theme_classic() +
    ylab(ylab) +
    xlab("") +
    theme(legend.position = "none")
}

#### WY ####
## load in composition data, clean and modify columns as usual
comp.wy <- read.csv("data/hpg_total.csv") #Wyoming species comp data 
comp.wy <- comp.wy %>% filter(year != "2020") #keep only 2023 data 
comp.wy$subplot <- as.factor(comp.wy$subplot)
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy <- comp.wy %>% unite(plot, c(block, trt, subplot), sep = ".", remove=F) # make unique plot variable
fornativecover <- comp.wy %>% filter(species!="BG"&
                                       species!="Litter"&
                                       native == "N") %>% #only native live veg
  group_by(year,block,trt,subplot) %>% 
  summarize(nativecov = sum(cover, na.rm=T)) #summarize total live veg per subplot
comp.wy <- merge(comp.wy,fornativecover, all.x = T)
comp.wy.wide <- comp.wy %>% select(-c("prob","native","graminoid")) #columns to drop 
comp.wy.wide <- comp.wy.wide %>% pivot_wider(id_cols = c("year","block","trt","subplot","drought",
                                                         "nativecov","BG", "Litter","plot","sub.tveg"), 
                                             names_from = "species", 
                                             values_from = "cover")
#comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data
#clean environment
rm(fornativecover)

###how many WY 2022+2023 data need to be dropped from CWM calculations (until we get trait data)
subwy <- comp.wy.wide[,-c(18:66)] %>% 
  mutate(propnative = nativecov/sub.tveg*100) %>% 
  filter(propnative < 80)
table(subwy$year)
table(comp.wy.wide$year)
181/512*100
#clean environment
rm(subwy)
rm(comp.wy.wide)

# must load newSelectSpecies function
# source("code/daniels_code/newSelectSpecies.R")

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

### pca
pca <- princomp(traits.wy[,3:9], cor=TRUE)
summary(pca)
biplot(pca)
traits.wy$pc1 <- pca$scores[,1]
traits.wy$pc2 <- pca$scores[,2]

#set color scheme
colors = brewer.pal(12,"Paired")
traits.wy$cols = c(colorRampPalette(colors)(nrow(traits.wy)))

### Subset monocots and dicots
grams <- subset(traits.wy, graminoid==1)
forbs <- subset(traits.wy, graminoid==0)

# Arrange trait matrix alphabetically
trait.matrix.wy <- traits.wy[order(rownames(traits.wy)),]
#Select colons of interest
trait.matrix.wy <- as.matrix(trait.matrix.wy[,2:11])

## pretreatment/ seeding probability 
# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
preds.wy <- read.csv("data/allplot.assemblages.csv") #data
preds.wy <- preds.wy %>% arrange(trt,block)
comms_p.wy <- labdsv::matrify(data.frame(preds.wy$trt,preds.wy$species,preds.wy$prob))
comms_p.wy <- comms_p.wy[,order(colnames(comms_p.wy))]
cwm_p.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms_p.wy), bin.num=c("graminoid","veg","c4"))

## 2021
#Arrange cover estimates for field data by year
hpg21 <- comp.wy %>% filter(year == "2021") %>% arrange(trt,block)%>%
  filter(native=="N") %>% filter(species != "CADU")
hpg21 <- hpg21 %>% unite(trt.b.sub, c(trt, block, subplot), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
comms21.wy <- labdsv::matrify(data.frame(hpg21$trt.b.sub,hpg21$species,hpg21$cover))
comms21.wy <- comms21.wy[,order(colnames(comms21.wy))]

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm21.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms21.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64*2
groups.wy <- c(rep("dt",N), rep("fd",N), rep("ir",N), rep("rand",N))

#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms21.wy <- vegan::metaMDS(cwm21.wy, distance="euclidean")
nms21.wy
plot(nms21.wy, type="t", main="euclidean")
ordiellipse(nms21.wy, groups.wy, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms21.wy,cwm21.wy)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm21.wy$trt <- factor(groups.wy)
# cwm_p$trt <- factor(groups)
cwm21.wy$trt <- factor(cwm21.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm21.wy$Treatments <- cwm21.wy$trt
levels(cwm21.wy$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#Calculating functional diversity of each trait and extracting RaoQ on L115
raoq <- FD::dbFD(as.matrix(trait.matrix.wy), as.matrix(comms21.wy))
leafn <- FD::dbFD(as.matrix(trait.matrix.wy[,"leafn"]), as.matrix(comms21.wy))
lop <- FD::dbFD(as.matrix(trait.matrix.wy[,"lop"]), as.matrix(comms21.wy))
ldmc <- FD::dbFD(as.matrix(trait.matrix.wy[,"ldmc"]), as.matrix(comms21.wy))
srl <- FD::dbFD(as.matrix(trait.matrix.wy[,"srl"]), as.matrix(comms21.wy))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.wy[,"rootdiam"]), as.matrix(comms21.wy))
veg <- FD::dbFD(as.matrix(trait.matrix.wy[,"veg"]), as.matrix(comms21.wy))
cwm_roaq21.wy <- data.frame(leafn=leafn$RaoQ,
                       lop=lop$RaoQ,
                       ldmc=ldmc$RaoQ,
                       srl=srl$RaoQ,
                       rootdiam=rootdiam$RaoQ,
                       veg=veg$RaoQ,
                       full=raoq$RaoQ,
                       trt=factor(groups.wy))
cwm_roaq21.wy$trt <- factor(cwm_roaq21.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm21.wy$block <- as.factor(matrix(unlist(strsplit(rownames(cwm21.wy)," ")),ncol=2, byrow=T)[,2])
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm21.wy <- cwm21.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                              !block %in% covered ~ "cntl")) 
## figures
leafn.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$leafn,"Leaf N") +
  geom_point(aes(y=quantile(cwm21.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).
test <- CWM_trait_FD_plot(cwm_roaq21.wy, cwm_roaq21.wy$leafn, "Leaf N")  

srl.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$srl,"Specific root length") +
  geom_point(aes(y=quantile(cwm21.wy$srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5)
  
veg.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$veg,"Vegetative spread potential") +
  ylim(0,1)

ldmc.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(cwm21.wy$ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) 
  
lop.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(cwm21.wy$lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5)

rootdiam.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$rootdiam,"Root diameter")

FD.wy.21 <- cwm_roaq21.wy %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,max(cwm_roaq21.wy$full+3)) +
  theme(legend.position  = "none")


## 2022
#Arrange cover estimates for field data by year
hpg22 <- comp.wy %>% filter(year == "2022") %>% arrange(trt,block) %>%
  filter(native=="N") %>% filter(species != "CADU")
  #filter(species %in% moretraits$species)
hpg22 <- hpg22 %>% unite(trt.b.sub, c(trt, block, subplot), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
comms22.wy <- labdsv::matrify(data.frame(hpg22$trt.b.sub,hpg22$species,hpg22$cover))
comms22.wy <- comms22.wy[,order(colnames(comms22.wy))]

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm22.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms22.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64*2
groups.wy <- c(rep("dt",N), rep("fd",N), rep("ir",N), rep("rand",N))

##plots with all 0 need to be removed to run code, but this will change when we get weed trait data
cwm22.wy <- cwm22.wy[rowSums(is.na(cwm22.wy)) != ncol(cwm22.wy), ]
#right now, changing cwm data should be fine, but if these data are needed later, use removed22.wy to progogate figures. 
#removed22 <- cwm22.wy[rowSums(is.na(cwm22.wy)) == ncol(cwm22.wy), ] # 3fd, 7ir, 4r (see removed plots)
groups.wy <-c(rep("dt",N), rep("fd",N-3), rep("ir",N-7), rep("rand",N-4)) # groups removing uncalculatable rows
  
#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms22.wy <- vegan::metaMDS(cwm22.wy, distance="euclidean")
nms22.wy
plot(nms22.wy, type="t", main="euclidean")
ordiellipse(nms22.wy, groups.wy, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms22.wy,cwm22.wy)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm22.wy$trt <- factor(groups.wy)
# cwm_p$trt <- factor(groups)
cwm22.wy$trt <- factor(cwm22.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm22.wy$Treatments <- cwm22.wy$trt
levels(cwm22.wy$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#remove communities that do not have any species data in 2022
comms22.wy <- comms22.wy %>%
  filter(rowSums(comms22.wy) != 0)

#Calculating functional diversity of each trait and extracting RaoQ on L115
raoq <- FD::dbFD(as.matrix(trait.matrix.wy), as.matrix(comms22.wy))
leafn <- FD::dbFD(as.matrix(trait.matrix.wy[,"leafn"]), as.matrix(comms22.wy))
lop <- FD::dbFD(as.matrix(trait.matrix.wy[,"lop"]), as.matrix(comms22.wy))
ldmc <- FD::dbFD(as.matrix(trait.matrix.wy[,"ldmc"]), as.matrix(comms22.wy))
srl <- FD::dbFD(as.matrix(trait.matrix.wy[,"srl"]), as.matrix(comms22.wy))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.wy[,"rootdiam"]), as.matrix(comms22.wy))
veg <- FD::dbFD(as.matrix(trait.matrix.wy[,"veg"]), as.matrix(comms22.wy))
cwm_roaq22.wy <- data.frame(leafn=leafn$RaoQ,
                       lop=lop$RaoQ,
                       ldmc=ldmc$RaoQ,
                       srl=srl$RaoQ,
                       rootdiam=rootdiam$RaoQ,
                       veg=veg$RaoQ,
                       full=raoq$RaoQ,
                       trt=factor(groups.wy))
cwm_roaq22.wy$trt <- factor(cwm_roaq22.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm22.wy$block <- as.factor(matrix(unlist(strsplit(rownames(cwm22.wy)," ")),ncol=2, byrow=T)[,2])
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm22.wy <- cwm22.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                              !block %in% covered ~ "cntl")) 

## figures
leafn.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$leafn,"Leaf N") +
  geom_point(aes(y=quantile(cwm22.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

srl.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$srl,"Specific root length") +
  geom_point(aes(y=quantile(cwm22.wy$srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5)

veg.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$veg,"Vegetative spread potential") +
  ylim(0,1)

ldmc.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(cwm22.wy$ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) 

lop.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(cwm22.wy$lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5)

rootdiam.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$rootdiam,"Root diameter")

FD.wy.22 <- cwm_roaq22.wy %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,max(cwm_roaq22.wy$full+3)) +
  theme(legend.position  = "none")

## 2023
#Arrange cover estimates for field data by year
hpg23 <- comp.wy %>% filter(year == "2023") %>% arrange(trt,block) %>%
  filter(native=="N") %>% filter(species != "CADU")
hpg23 <- hpg23 %>% unite(trt.b.sub, c(trt, block, subplot), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
comms23.wy <- labdsv::matrify(data.frame(hpg23$trt.b.sub,hpg23$species,hpg23$cover))
comms23.wy <- comms23.wy[,order(colnames(comms23.wy))]

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm23.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms23.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64*2
groups <- c(rep("dt",N), rep("fd",N), rep("ir",N), rep("rand",N))

## REMOVE
# ##plots with all 0 need to be removed to run code, but this will change when we get weed trait data
# cwm23.wy <- cwm23.wy[rowSums(is.na(cwm23.wy)) != ncol(cwm23.wy), ] 
# #right now, changing cwm data should be fine, but if these data are needed later, use removed23.wy to propagate figures. 
# #removed23.wy <- cwm23.wy[rowSums(is.na(cwm23.wy)) == ncol(cwm23.wy), ] # 3fd, 7ir, 4r (see removed plots)
# groups.wy <-c(rep("dt",N), rep("fd",N-3), rep("ir",N-7), rep("rand",N-4)) # groups removing uncalculatable rows

#nonmetric multidimensional scaling and ploting of ellipses by treatment
nms23.wy <- vegan::metaMDS(cwm23.wy, distance="euclidean")
nms23.wy
plot(nms23.wy, type="t", main="euclidean")
ordiellipse(nms23.wy, groups, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms23.wy,cwm23.wy)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm23.wy$trt <- factor(groups)
# cwm_p$trt <- factor(groups)
cwm23.wy$trt <- factor(cwm23.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
cwm23.wy$Treatments <- cwm23.wy$trt
levels(cwm23.wy$Treatments) <- c("Invasion Resistant","Drought Tolerant","Functionally Diverse","Random")

#remove species not occurring in any community this year from comms and trait matrix
colSums(comms23.wy) #find species
comms23.wy <- comms23.wy %>% select(-ARPU)
remove23.wy <- "ARPU"
trait.matrix.wy21 <- trait.matrix.wy[!row.names(trait.matrix.wy)%in%remove23.wy,]

#Calculating functional diversity of each trait and extracting RaoQ on L115
raoq <- FD::dbFD(as.matrix(trait.matrix.wy21), as.matrix(comms23.wy))
leafn <- FD::dbFD(as.matrix(trait.matrix.wy21[,"leafn"]), as.matrix(comms23.wy))
lop <- FD::dbFD(as.matrix(trait.matrix.wy21[,"lop"]), as.matrix(comms23.wy))
ldmc <- FD::dbFD(as.matrix(trait.matrix.wy21[,"ldmc"]), as.matrix(comms23.wy))
srl <- FD::dbFD(as.matrix(trait.matrix.wy21[,"srl"]), as.matrix(comms23.wy))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.wy21[,"rootdiam"]), as.matrix(comms23.wy))
veg <- FD::dbFD(as.matrix(trait.matrix.wy21[,"veg"]), as.matrix(comms23.wy))
cwm_roaq23.wy <- data.frame(leafn=leafn$RaoQ,
                       lop=lop$RaoQ,
                       ldmc=ldmc$RaoQ,
                       srl=srl$RaoQ,
                       rootdiam=rootdiam$RaoQ,
                       veg=veg$RaoQ,
                       full=raoq$RaoQ,
                       trt=factor(groups))
cwm_roaq23.wy$trt <- factor(cwm_roaq23.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm23.wy$block <- as.factor(matrix(unlist(strsplit(rownames(cwm23.wy)," ")),ncol=2, byrow=T)[,2])
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm23.wy <- cwm23.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                              !block %in% covered ~ "cntl")) 

## figures
leafn.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$leafn,"Leaf N") +
  geom_point(aes(y=quantile(cwm23.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

srl.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$srl,"Specific root length") +
  geom_point(aes(y=quantile(cwm23.wy$srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5)

veg.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$veg,"Vegetative spread potential") +
  ylim(0,1)

ldmc.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(cwm23.wy$ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) 

lop.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(cwm23.wy$lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5)

rootdiam.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$rootdiam,"Root diameter")

FD.wy.23 <- cwm_roaq23.wy %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,max(cwm_roaq23.wy$full+3)) +
  theme(legend.position  = "none")

#FD of traits if desired:
#leafn.wy <- CWM_trait_plot(cwm21.wy,leafn,traits.wy)

## merge cwm's together for storing
wpre <- cwm21.wy %>% mutate(year = "0")
w21 <- cwm21.wy %>% mutate(year = "2021") #add year
w22 <- cwm21.wy %>% mutate(year = "2022") #add year
w23 <- cwm21.wy %>% mutate(year = "2023") #add year

cwm.wy <- bind_rows(wpre, w21) #bind 1st
cwm.wy <- bind_rows(cwm.wy, w22) #bind again
cwm.wy <- bind_rows(cwm.wy, w23) #bind again

#save 
write.csv(cwm.wy, "data/cwm_wy.csv", row.names = F)

## Produce a .tiff file of our plots together!
library(patchwork)
# export 21
tiff("figures/cwm WY/traits_2021.tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.wy.21 + srl.wy.21 + veg.wy.21) / (ldmc.wy.21 + lop.wy.21 + rootdiam.wy.21)
dev.off()
tiff("figures/cwm WY/FD_2021.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.wy.21
dev.off()
# export 22
tiff("figures/cwm WY/traits_2022.tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.wy.22 + srl.wy.22 + veg.wy.22) / (ldmc.wy.22 + lop.wy.22 + rootdiam.wy.22)
dev.off()
tiff("figures/cwm WY/FD_2022.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.wy.22
dev.off()
# export 23
tiff("figures/cwm WY/traits_2023.tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.wy.23 + srl.wy.23 + veg.wy.23) / (ldmc.wy.23 + lop.wy.23 + rootdiam.wy.23)
dev.off()
tiff("figures/cwm WY/FD_2023.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.wy.23
dev.off()

# tiff("figures/cwm/traits_2020.tiff", width=8.5, height=6, "in", res=400, compression = "lzw")
# gridExtra::grid.arrange(a,b,c,d,e,f,ncol=3)
# dev.off()
# tiff("figures/trait_raoq_2023.tiff",  width=8.5, height=6, "in", res=400, compression = "lzw")
# gridExtra::grid.arrange(a,b,c,d,e,f,ncol=3)
# dev.off()
# 
# tiff("figures/multivariate_raoq_2023.tiff",  width=4, height=4, "in", res=400, compression = "lzw")
# cwm_roaq %>% 
#   ggplot(aes(trt,full)) +
#   geom_violin(aes(fill=trt), width=1, trim=F) +
#   geom_boxplot(width=.25) +
#   geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
#   scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
#   theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
#   ylim(0,max(cwm_roaq$full+3)) +
#   theme(legend.position  = "none")
# 
# dev.off()


#### CA ####
## load in composition data, clean and modify columns as usual
comp.ca <- read.csv("data/Species_Composition_allyears.csv") #read in California comp data
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
# remove species not present in other comp dataframe
trait.matrix.ca <- trait.matrix.ca %>% 
  filter(Code != "AVBA") %>% 
  filter(Code != "BRNI") %>%
  filter(Code != "CAME") %>%
  #filter(Code != "FEMY") %>% 
  filter(Code != "LUAL") %>%
  filter(Code != "MASA") %>% 
  filter(Code != "PEHE") %>% 
  filter(Code != "SACO") 

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
comms_p.ca <- comms_p.ca %>% select(-CAME)
trait.matrix.ca.pred <- traits.ca[order(rownames(traits.ca)),]
# remove species not present in other comp dataframe
trait.matrix.ca.pred <- trait.matrix.ca.pred %>% 
  filter(Code != "AVBA") %>% 
  filter(Code != "BRNI") %>%
  filter(Code != "LUAL") %>%
  filter(Code != "MASA") %>% 
  filter(Code != "PEHE") %>% 
  filter(Code != "SACO")
test <- data.frame(trait.matrix.ca.pred[,2:8], row.names = trait.matrix.ca.pred[,1])
rownames(test)
trait.matrix.ca.pred <- as.matrix(test)
trait.matrix.ca.pred <- trait.matrix.ca.pred[order(rownames(trait.matrix.ca.pred)),]
remove22.pred <- c("FEMY","FEPE","BRMA","CAME")
trait.matrix.ca.pred <- trait.matrix.ca.pred[!row.names(trait.matrix.ca.pred)%in%remove22.pred,]
cwm_p.ca <- FD::functcomp(as.matrix(trait.matrix.ca.pred), as.matrix(comms_p.ca), bin.num=c("graminoid"))

## pretreatment (not done)
#Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
# preds <- preds %>% arrange(trt,block)
# comms_p <- labdsv::matrify(data.frame(preds$trt,preds$species,preds$prob))
# comms_p <- comms_p[,order(colnames(comms_p))]
# cwm_p <- FD::functcomp(as.matrix(trait.matrix.ca), as.matrix(comms_p), bin.num=c("graminoid","veg","c4"))

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

# remove columns from comms (communties) that are not present in trait.matrix (can we get these? too small proportion?)
comms21.ca <- comms21.ca %>% select(-c("CACI", "BRHO", "BRDI", "HOMU"))
#remove species from matrix with NA in all communities this year (not present yet) 
trait.matrix.ca21 <- trait.matrix.ca[!row.names(trait.matrix.ca)=="FEMY",]
# rownames(trait.matrix.ca22)
# colnames(comms21.ca)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm21.ca <- FD::functcomp(as.matrix(trait.matrix.ca21), as.matrix(comms21.ca), bin.num=c("graminoid"))

# Building out treatment identification which was absent from our original csv.
N <- 53
groups.ca <- c(rep("dt",N), rep("fd",N-1), rep("ir",N-1), rep("rand",N))

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
trait.matrix.ca21 <- trait.matrix.ca21[!row.names(trait.matrix.ca21)%in%remove21,]

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
cwm21.ca$block <- as.factor(matrix(unlist(strsplit(rownames(cwm21.ca),"[.]")),ncol=2, byrow=T)[,2])
# Define drought treatment at block level
block.water <- hpf21 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cwm21.ca <- merge(cwm21.ca, block.water, by = c("block","trt"), all.x=T) #merge


## figures
leafn.ca.21 <- CWM_trait_plot(cwm21.ca,cwm21.ca$N,"Leaf N") +
  geom_point(aes(y=quantile(cwm21.ca$N,.25),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

srl.ca.21 <- CWM_trait_plot(cwm21.ca,cwm21.ca$SRL,"Specific root length") +
  geom_point(aes(y=quantile(cwm21.ca$SRL,.7557),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5)

rmf.ca.21 <- CWM_trait_plot(cwm21.ca,cwm21.ca$RMF,"Root mass fraction") +
  geom_point(aes(y=quantile(cwm21.ca$RMF,.25),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

lma.ca.21 <- CWM_trait_plot(cwm21.ca,cwm21.ca$LMA,"Leaf mass per area") +
  geom_point(aes(y=quantile(cwm21.ca$LMA,.75),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) 

seedmass.ca.21 <- CWM_trait_plot(cwm21.ca,cwm21.ca$seed.mass,"Seed mass") +
  geom_point(aes(y=quantile(cwm21.ca$seed.mass,.75),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5)

rootdiam.ca.21 <- CWM_trait_plot(cwm21.ca,cwm21.ca$Rdiam,"Root diameter") #FD

FD.ca.21 <- cwm_roaq21.ca %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,max(cwm_roaq21.ca$full+3)) +
  theme(legend.position  = "none")


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

# remove columns from comms (communties) that are not present in trait.matrix (can we get these? too small proportion?)
comms22.ca <- comms22.ca %>% select(-c("CACI", "BRDI"))
#remove species from matrix with NA in all communities this year (not present yet) 
trait.matrix.ca22 <- trait.matrix.ca[!row.names(trait.matrix.ca)=="FEMY",]
 rownames(trait.matrix.ca22)
 colnames(comms22.ca)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm22.ca <- FD::functcomp(as.matrix(trait.matrix.ca22), as.matrix(comms22.ca), bin.num=c("graminoid"))

# Building out treatment identification which was absent from our original csv.
N <- 53
# plots with all 0 (no spp in community) need to be removed to run nms2,
cwm22.ca <- cwm22.ca[rowSums(is.na(cwm22.ca)) != ncol(cwm22.ca), ] 
#right now, changing cwm22.ca data should be fine, but if these data are needed later, use removed23.wy to propagate figures. 
#removed22.ca <- cwm22.ca[rowSums(is.na(cwm22.ca)) == ncol(cwm22.ca), ] # fd.36 and ir.47 (see removed plots)
groups.ca <-c(rep("dt",N), rep("fd",N-2), rep("ir",N-2), rep("rand",N)) # groups removing uncalculatable rows

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
trait.matrix.ca22 <- trait.matrix.ca22[!row.names(trait.matrix.ca22)%in%remove22.ca,]

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
cwm22.ca$block <- as.factor(matrix(unlist(strsplit(rownames(cwm22.ca),"[.]")),ncol=2, byrow=T)[,2])
# Define drought treatment at block level
block.water <- hpf22 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cwm22.ca <- merge(cwm22.ca, block.water, by = c("block","trt"), all.x=T) #merge


## figures
leafn.ca.22 <- CWM_trait_plot(cwm22.ca,cwm22.ca$N,"Leaf N") +
  geom_point(aes(y=quantile(cwm22.ca$N,.25),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

srl.ca.22 <- CWM_trait_plot(cwm22.ca,cwm22.ca$SRL,"Specific root length") +
  geom_point(aes(y=quantile(cwm22.ca$SRL,.7557),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5)

rmf.ca.22 <- CWM_trait_plot(cwm22.ca,cwm22.ca$RMF,"Root mass fraction") +
  geom_point(aes(y=quantile(cwm22.ca$RMF,.25),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

lma.ca.22 <- CWM_trait_plot(cwm22.ca,cwm22.ca$LMA,"Leaf mass per area") +
  geom_point(aes(y=quantile(cwm22.ca$LMA,.75),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) 

seedmass.ca.22 <- CWM_trait_plot(cwm22.ca,cwm22.ca$seed.mass,"Seed mass") +
  geom_point(aes(y=quantile(cwm22.ca$seed.mass,.75),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5)

rootdiam.ca.22 <- CWM_trait_plot(cwm22.ca,cwm22.ca$Rdiam,"Root diameter") #FD

FD.ca.22 <- cwm_roaq22.ca %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,max(cwm_roaq22.ca$full+3)) +
  theme(legend.position  = "none")


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

# remove columns not present in communities this year
comms23.ca <- comms23.ca %>% select(-c("CACI", "BRDI", "BRHO", "HOMU", "POMO"))
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
groups.ca <-c(rep("dt",N), rep("fd",N-2), rep("ir",N-1), rep("rand",N)) # groups removing uncalculatable rows

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
cwm23.ca$block <- as.factor(matrix(unlist(strsplit(rownames(cwm23.ca),"[.]")),ncol=2, byrow=T)[,2])
# Define drought treatment at block level
block.water <- hpf23 %>% select(c(block,trt,water))# get data from OG dataframe
block.water$trt <- tolower(block.water$trt) #make lowercase to match
block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
cwm23.ca <- merge(cwm23.ca, block.water, by = c("block","trt"), all.x=T) #merge


## figures
leafn.ca.23 <- CWM_trait_plot(cwm23.ca,cwm23.ca$N,"Leaf N") +
  geom_point(aes(y=quantile(cwm23.ca$N,.25),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

srl.ca.23 <- CWM_trait_plot(cwm23.ca,cwm23.ca$SRL,"Specific root length") +
  geom_point(aes(y=quantile(cwm23.ca$SRL,.7557),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5)

rmf.ca.23 <- CWM_trait_plot(cwm23.ca,cwm23.ca$RMF,"Root mass fraction") +
  geom_point(aes(y=quantile(cwm23.ca$RMF,.25),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) # Specifying the target object (red dot).

lma.ca.23 <- CWM_trait_plot(cwm23.ca,cwm23.ca$LMA,"Leaf mass per area") +
  geom_point(aes(y=quantile(cwm23.ca$LMA,.75),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) 

seedmass.ca.23 <- CWM_trait_plot(cwm23.ca,cwm23.ca$seed.mass,"Seed mass") +
  geom_point(aes(y=quantile(cwm23.ca$seed.mass,.75),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5)

rootdiam.ca.23 <- CWM_trait_plot(cwm23.ca,cwm23.ca$Rdiam,"Root diameter") #FD

FD.ca.23 <- cwm_roaq23.ca %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,max(cwm_roaq23.ca$full+3)) +
  theme(legend.position  = "none")

#FD of traits if desired:
#leafn.wy <- CWM_trait_plot(cwm21.wy,leafn,traits.wy)

## merge cwm's together for storing
cpre <- cwm21.ca %>% mutate(year = "0") #add prob
c21 <- cwm21.ca %>% mutate(year = "2021") #add year
c22 <- cwm21.ca %>% mutate(year = "2022") #add year
c23 <- cwm21.ca %>% mutate(year = "2023") #add year

cwm.ca <- bind_rows(cpre, c21) #bind 1st
cwm.ca <- bind_rows(cwm.ca, c22) #bind again
cwm.ca <- bind_rows(cwm.ca, c23) #bind again

#save 
write.csv(cwm.ca, "data/cwm_ca.csv", row.names = F)

## Produce a .tiff file of our plots together!
library(patchwork)
# export 21
tiff("figures/cwm CA/traits_2021.tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.ca.21 + srl.ca.21 + rmf.ca.21) / (lma.ca.21 + seedmass.ca.21 + rootdiam.ca.21)
dev.off()
tiff("figures/cwm CA/FD_2021.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.ca.21
dev.off()
# export 22
tiff("figures/cwm CA/traits_2022.tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.ca.22 + srl.ca.22 + rmf.ca.22) / (lma.ca.22 + seedmass.ca.22 + rootdiam.ca.22)
dev.off()
tiff("figures/cwm CA/FD_2022.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.ca.22
dev.off()
# export 23
tiff("figures/cwm CA/traits_2023.tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.ca.23 + srl.ca.23 + rmf.ca.23) / (lma.ca.23 + seedmass.ca.23 + rootdiam.ca.23)
dev.off()
tiff("figures/cwm CA/FD_2023.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.ca.23
dev.off()


#### add a legent to plots?
# leg <- cwm %>% 
#   ggplot(aes(Treatments,rootdiam)) +
#   geom_violin(aes(fill=Treatments), width=1, trim=F) +
#   scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
#   theme_classic()
# leg
# 
# tiff("trait_legend.tiff", width = 3, height = 2, "in", res=400, compression = "lzw")
# grid::grid.newpage()
# grid::grid.draw(cowplot::get_legend(leg))
# dev.off()
# 
# gridExtra::grid.arrange(a,b,c,d,e,f,ncol=3)