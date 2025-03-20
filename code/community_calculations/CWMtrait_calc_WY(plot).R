#### Calculate the CWM traits for each plot in each year based on the 
#### taxonomic composition for Wyoming site. 

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


#### WY ####
## load in composition data, clean and modify columns as usual
comp.wy <- read.csv("data/comp_wy_plot.csv") #Wyoming species comp data PLOT level
comp.wy <- comp.wy %>% filter(year != "2020") #keep only 2023 data 
comp.wy$drought <- as.factor(comp.wy$drought)
comp.wy$year<-as.factor(comp.wy$year)
comp.wy$trt <- as.factor(comp.wy$trt)
comp.wy$block <- as.factor(comp.wy$block)
comp.wy <- comp.wy %>% unite(plot, c(block, trt), sep = ".", remove=F) # make unique plot variable

###how many WY 2022+2023 data need to be dropped from CWM calculations
subwy <- comp.wy %>%
  mutate(propnative = nativecov.plot/totcov.plot*100) %>%
  filter(propnative < 80)
table(subwy$year)
table(comp.wy$year)
(17+132)/(256*3)*100 #losing 19% total across years

# must load newSelectSpecies function
# source("code/daniels_code/newSelectSpecies.R")

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
og25 <- unique(colnames(comms_p.wy)) #make list of 25 native species
cwm_p.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms_p.wy), bin.num=c("graminoid","veg","c4"))
#Define block by extracting the numeric from the cwm rownames
cwm_p.wy$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(cwm_p.wy)))
#Define trt by extracting the subplot from the cwm rownames
cwm_p.wy$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(cwm_p.wy))))
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm_p.wy <- cwm_p.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                                    !block %in% covered ~ "cntl")) 

#Calculating functional diversity of each trait and extracting RaoQ (RoaQ is not scaled here, needs to be done when using these columns)
raoq <- FD::dbFD(as.matrix(trait.matrix.wy), as.matrix(comms_p.wy))
rootdiam <- FD::dbFD(as.matrix(trait.matrix.wy[,"rootdiam"]), as.matrix(comms_p.wy))
veg <- FD::dbFD(as.matrix(trait.matrix.wy[,"veg"]), as.matrix(comms_p.wy))
cwm_roaq_p.wy <- data.frame(rootdiam=rootdiam$RaoQ,
                            veg=veg$RaoQ,
                            full=raoq$RaoQ)
#cwm_roaq_p.wy$trt <- factor(cwm_roaq_p.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)


## 2021
#Arrange cover estimates for field data by year and remove weed species for now
hpg21 <- comp.wy %>% filter(year == "2021") %>% arrange(trt,block)
hpg21 <- hpg21 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
needcols <- c(og25,"trt.b")
comms21.wy <- hpg21 %>% column_to_rownames("trt.b")
comms21.wy <- comms21.wy %>% select(all_of(og25)) #ONLY NATIVE 25 and grouping columns
comms21.wy <- comms21.wy[,order(colnames(comms21.wy))]
comms21.wy <- replace(comms21.wy, is.na(comms21.wy), 0)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm21.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms21.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64*1
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

#Calculating functional diversity of each trait and extracting RaoQ (RoaQ is not scaled here, needs to be done when using these columns)
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
cwm21.wy$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm21.wy)))
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm21.wy <- cwm21.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                                    !block %in% covered ~ "cntl")) 

## 2022
#Arrange cover estimates for field data by year and remove weed species for now
hpg22 <- comp.wy %>% filter(year == "2022") %>% arrange(trt,block)
hpg22 <- hpg22 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
comms22.wy <- hpg22 %>% column_to_rownames("trt.b")
comms22.wy <- comms22.wy %>% select(all_of(og25)) #ONLY NATIVE 25 and grouping columns
comms22.wy <- comms22.wy[,order(colnames(comms22.wy))]
comms22.wy <- replace(comms22.wy, is.na(comms22.wy), 0)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm22.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms22.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64
#groups.wy <- c(rep("dt",N), rep("fd",N), rep("ir",N), rep("rand",N))

##plots with all 0 need to be removed to run code, but this will change when we get weed trait data
cwm22.wy <- cwm22.wy[rowSums(is.na(cwm22.wy)) != ncol(cwm22.wy), ]
#right now, changing cwm data should be fine, but if these data are needed later, use removed22.wy to progogate figures. 
#removed22 <- cwm22.wy[rowSums(is.na(cwm22.wy)) == ncol(cwm22.wy), ] # ir12 and ir18 (see removed plots)
groups.wy <-c(rep("dt",N), rep("fd",N), rep("ir",N-2), rep("rand",N)) # groups removing uncalculatable rows

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
cwm22.wy$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm22.wy)))
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm22.wy <- cwm22.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                                    !block %in% covered ~ "cntl")) 


## 2023
#Arrange cover estimates for field data by year and remove weed species for now
hpg23 <- comp.wy %>% filter(year == "2023") %>% arrange(trt,block)
hpg23 <- hpg23 %>% unite(trt.b, c(trt, block), sep = ".", remove=F) # make unique plot variable

# Reshape field data into a matrix and arrange alphabetically
comms23.wy <- hpg23 %>% column_to_rownames("trt.b")
comms23.wy <- comms23.wy %>% select(all_of(og25)) #ONLY NATIVE 25 and grouping columns
comms23.wy <- comms23.wy[,order(colnames(comms23.wy))]
comms23.wy <- replace(comms23.wy, is.na(comms23.wy), 0)

#Calculate community weighted means. Note that bin.num must bne specified for binary outcomes
cwm23.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms23.wy), bin.num=c("graminoid","veg","c4"))

# Building out treatment identification which was absent from our original csv.
N <- 64
groups.wy <- c(rep("dt",N), rep("fd",N), rep("ir",N), rep("rand",N))

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
ordiellipse(nms23.wy, groups.wy, col=c("orange","green","red","darkgray"), conf=0.95)
vectors <- envfit(nms23.wy,cwm23.wy)
plot(vectors,p.max=0.05)

# Factoring the groups and providing their full names for plotting
cwm23.wy$trt <- factor(groups.wy)
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
                            trt=factor(groups.wy))
cwm_roaq23.wy$trt <- factor(cwm_roaq23.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
#Define block by extracting the numeric from the cwm rownames
cwm23.wy$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(cwm23.wy)))
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm23.wy <- cwm23.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                                    !block %in% covered ~ "cntl")) 


## merge cwm's together for storing
wpre <- cwm_p.wy %>% mutate(year = "0") 
wpre$trt <- as.factor(as.character(wpre$trt))
wpre <- wpre %>% mutate(trt = str_replace(trt, "^r $", "rand")) #make r match rand 
wpre <- wpre %>% mutate(trt = str_replace(trt, "^ir $", "ir")) #make r match rand 
wpre <- wpre %>% mutate(trt = str_replace(trt, "^fd $", "fd")) #make r match rand 
wpre <- wpre %>% mutate(trt = str_replace(trt, "^dt $", "dt")) #make r match rand 
#wpre$trt <- factor(wpre$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE)
rownames(wpre) <- NULL
w21 <- cwm21.wy %>% mutate(year = "2021") #add year
w21$trt <- as.factor(as.character(w21$trt))
rownames(w21) <- NULL
w22 <- cwm22.wy %>% mutate(year = "2022") #add year
w22$trt <- as.factor(as.character(w22$trt))
rownames(w22) <- NULL
w23 <- cwm23.wy %>% mutate(year = "2023") #add year
w23$trt <- as.factor(as.character(w23$trt))
rownames(w23) <- NULL

cwm.wy <- bind_rows(wpre, w21) #bind 1st
cwm.wy <- bind_rows(cwm.wy, w22) #bind again
cwm.wy <- bind_rows(cwm.wy, w23) #bind again
#save 
write.csv(cwm.wy, "data/cwm_wy(plot).csv", row.names = F)

## merge cwm's together for storing
## I am only calculating RaoQ for preds for traits I need it, this is the only reason NA appears in the dataframe
wfd0 <- cwm_roaq_p.wy %>% mutate(year = "0") #add year
wfd0$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(wfd0))) #Define block by extracting the numeric from the cwm rownames
wfd0$trt <- as.factor(sub("([a-z]+)\\s(\\d+)$", "\\1", rownames(wfd0))) #Define trt by extracting the trt from the cwm rownames
wfd0 <- wfd0 %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand 
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
wfd0 <- wfd0 %>% mutate(drought = case_when(block %in% covered ~ "drt", # Define drought treatment at block level
                                            !block %in% covered ~ "cntl"))
wfd1 <- cwm_roaq21.wy %>% mutate(year = "2021") #add year
wfd1$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(wfd1))) #Define block by extracting the numeric from the cwm rownames
#wfd1$subplot <- as.factor(sub(".*\\.(\\d+)\\.(\\w)$", "\\2", rownames(wfd1))) #Define subplot by extracting the subplot from the cwm rownames
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
wfd1 <- wfd1 %>% mutate(drought = case_when(block %in% covered ~ "drt", # Define drought treatment at block level
                                            !block %in% covered ~ "cntl"))
wfd2 <- cwm_roaq22.wy %>% mutate(year = "2022") #add year
wfd2$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(wfd2))) #Define block by extracting the numeric from the cwm rownames
#wfd2$subplot <- as.factor(sub(".*\\.(\\d+)\\.(\\w)$", "\\2", rownames(wfd2))) #Define subplot by extracting the subplot from the cwm rownames
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
wfd2 <- wfd2 %>% mutate(drought = case_when(block %in% covered ~ "drt", # Define drought treatment at block level
                                            !block %in% covered ~ "cntl"))
wfd3 <- cwm_roaq23.wy %>% mutate(year = "2023") #add year
wfd3$block <- as.factor(sub(".*\\.(\\d+)$", "\\1", rownames(wfd3))) #Define block by extracting the numeric from the cwm rownames
#wfd3$subplot <- as.factor(sub(".*\\.(\\d+)\\.(\\w)$", "\\2", rownames(wfd3))) #Define subplot by extracting the subplot from the cwm rownames
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
wfd3 <- wfd3 %>% mutate(drought = case_when(block %in% covered ~ "drt", # Define drought treatment at block level
                                            !block %in% covered ~ "cntl"))


cwm.wyfd <- bind_rows(wfd0, wfd1) #bind 1st
cwm.wyfd <- bind_rows(cwm.wyfd, wfd2) #bind again
cwm.wyfd <- bind_rows(cwm.wyfd, wfd3) #bind again
#save 
write.csv(cwm.wyfd, "data/cwm_raoq_wy(plot).csv", row.names = F)



#### FIGURES ####
#### CWM ANOVA's and figures for WY site
#### How to CWM's differ by trt?
#### Creates plots to visualize mean, variance, and target for each traits:

#### Functions needed to normalizing RoaQ, making post-hoc letters, 
#### and creating figures
# Function for normalizing FD
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# automatically extract post-hoc comparisons from tukey and place on plot
generateTukeyLabel <- function(model, y_var) {
  # Perform Tukey HSD
  tukey_cov <- TukeyHSD(model)
  # Display letters for each site
  multcompView::multcompLetters4(model, tukey_cov)
  sum_labels <- data.frame(multcompView::multcompLetters4(model, tukey_cov)$trt['Letters'])
  sum_labels <- sum_labels %>% rownames_to_column("trt")
  # Obtain letter position for y axis using means
  yvalue <- aggregate(y_var ~ trt, data = model$model, max)
  # Merge dataframes
  tukeylabel <- merge(sum_labels, yvalue)
  return(tukeylabel)
}

# make violin plot for traits functional diversity (RoaQ)
CWM_trait_FD_plot <- function(data, trait, ylab) {
  ggplot(data, aes_string("trt", trait)) +
    geom_violin(aes(fill = trt), width = 1, trim = FALSE) +
    geom_boxplot(width = 0.25) +
    #geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.3) +
    scale_fill_viridis_d(option = "D", begin = 0.1, end = 1, alpha = 0.7) +
    theme_classic() +
    ylab(paste("FD (RaoQ) of ", ylab)) +
    xlab("") +
    ylim(0,1)+
    #ylim(0,max(trait+.5)) +
    theme(legend.position = "none")
}

# make violin plot for trait 
CWM_trait_plot <- function(data, trait, ylab) {
  ggplot(data, aes_string("trt", trait)) +
    geom_violin(aes(fill = trt), width = 1, trim = FALSE) +
    geom_boxplot(width = 0.25) +
    #geom_jitter(width = 0.12, height = 0, size = 1, alpha = 0.3) +
    scale_fill_viridis_d(option = "D", begin = 0.1, end = 1, alpha = 0.7) +
    theme_classic() +
    ylab(ylab) +
    xlab("") +
    theme(legend.position = "none")
}


#### reading in data
dat <- read.csv("data/cwm_wy(plot).csv") #cwm data
dat$trt <- factor(dat$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
datFD <- read.csv("data/cwm_raoq_wy(plot).csv") #cwm RoaQ data
datFD$trt <- factor(datFD$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
#trait data
traits.wy <- read.csv("data/trait_data/mixedgrass.csv", header=TRUE, row.names=1) # CSV of species-trait combinations (for OG 25)
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

## subset data by year for figures
subdat21 <- dat %>% filter(year=="2021")
subdatFD21 <- datFD %>% filter(year=="2021")
subdat22 <- dat %>% filter(year=="2022")
subdatFD22 <- datFD %>% filter(year=="2022")
subdat23 <- dat %>% filter(year=="2023")
subdatFD23 <- datFD %>% filter(year=="2023")

#### Did traits of seeding treatments differ? (ANOVA, post-hoc Tukey, and plotting)
#### 2021 ####
summary(leafn21.mod <- aov(leafn~trt, subdat21))
leafn21.tuk <- generateTukeyLabel(leafn21.mod, subdat21$leafn)
leafn.wy.21 <- CWM_trait_plot(subdat21,subdat21$leafn,"Leaf N") +
  geom_point(aes(y=quantile(leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
  #  geom_point(aes(y=quantile(cwm21.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl21.mod <- aov(srl~trt, subdat21))
srl21.tuk <- generateTukeyLabel(srl21.mod, subdat21$srl)
srl.wy.21 <- CWM_trait_plot(subdat21,subdat21$srl,"Specific root length") +
  geom_point(aes(y=quantile(srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = srl21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

# veg.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$veg,"Vegetative spread potential") +
#   ylim(0,1)
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD21$veg <- normalize(subdatFD21$veg) #(normalize FD/RoaQ)
summary(veg21.mod <- aov(srl~trt, subdatFD21)) #summary model
veg21.tuk <- generateTukeyLabel(veg21.mod, subdatFD21$veg) #tukey and labels
veg.wy.21 <- CWM_trait_FD_plot(subdatFD21,subdatFD21$veg,"Vegetative spread potential") +
  geom_point(aes(y=max(veg),x="ir"),data=subdatFD21, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = veg21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=2, col="black") #add tukey labels

summary(ldmc21.mod <- aov(ldmc~trt, subdat21))
ldmc21.tuk <- generateTukeyLabel(ldmc21.mod, subdat21$ldmc)
ldmc.wy.21 <- CWM_trait_plot(subdat21,subdat21$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = ldmc21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lop21.mod <- aov(lop~trt, subdat21))
lop21.tuk <- generateTukeyLabel(lop21.mod, subdat21$lop)
lop.wy.21 <- CWM_trait_plot(subdat21,subdat21$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = lop21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2, vjust=-1.75, col="black") #add tukey labels

# rootdiam.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD21$rootdiam <- normalize(subdatFD21$rootdiam) #(normalize FD/RoaQ)
summary(rd21.mod <- aov(rootdiam~trt, subdatFD21)) #summary model
rd21.tuk <- generateTukeyLabel(rd21.mod, subdatFD21$rootdiam) #tukey and labels
rootdiam.wy.21 <- CWM_trait_FD_plot(subdatFD21,subdatFD21$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=subdatFD21, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

subdatFD21$full <- normalize(subdatFD21$full) #(normalize FD/RoaQ)
summary(full21.mod <- aov(full~trt, subdatFD21)) #summary model
full21.tuk <- generateTukeyLabel(full21.mod, subdatFD21$full) #tukey and labels
FD.wy.21 <- subdatFD21 %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq21.wy$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=subdatFD21, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels



#### 2022 ####
summary(leafn22.mod <- aov(leafn~trt, subdat22))
leafn22.tuk <- generateTukeyLabel(leafn22.mod, subdat22$leafn)
leafn.wy.22 <- CWM_trait_plot(subdat22,subdat22$leafn,"Leaf N") +
  geom_point(aes(y=quantile(leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
  #  geom_point(aes(y=quantile(cwm22.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl22.mod <- aov(srl~trt, subdat22))
srl22.tuk <- generateTukeyLabel(srl22.mod, subdat22$srl)
srl.wy.22 <- CWM_trait_plot(subdat22,subdat22$srl,"Specific root length") +
  geom_point(aes(y=quantile(srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = srl22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

# veg.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$veg,"Vegetative spread potential") +
#   ylim(0,1)
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD22$veg <- normalize(subdatFD22$veg) #(normalize FD/RoaQ)
summary(veg22.mod <- aov(srl~trt, subdatFD22)) #summary model
veg22.tuk <- generateTukeyLabel(veg22.mod, subdatFD22$veg) #tukey and labels
veg.wy.22 <- CWM_trait_FD_plot(subdatFD22,subdatFD22$veg,"Vegetative spread potential") +
  geom_point(aes(y=max(veg),x="ir"),data=subdatFD22, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = veg22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=2, col="black") #add tukey labels

summary(ldmc22.mod <- aov(ldmc~trt, subdat22))
ldmc22.tuk <- generateTukeyLabel(ldmc22.mod, subdat22$ldmc)
ldmc.wy.22 <- CWM_trait_plot(subdat22,subdat22$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = ldmc22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lop22.mod <- aov(lop~trt, subdat22))
lop22.tuk <- generateTukeyLabel(lop22.mod, subdat22$lop)
lop.wy.22 <- CWM_trait_plot(subdat22,subdat22$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = lop22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2, vjust=-1.75, col="black") #add tukey labels

# rootdiam.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD22$rootdiam <- normalize(subdatFD22$rootdiam) #(normalize FD/RoaQ)
summary(rd22.mod <- aov(rootdiam~trt, subdatFD22)) #summary model
rd22.tuk <- generateTukeyLabel(rd22.mod, subdatFD22$rootdiam) #tukey and labels
rootdiam.wy.22 <- CWM_trait_FD_plot(subdatFD22,subdatFD22$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=subdatFD22, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

subdatFD22$full <- normalize(subdatFD22$full) #(normalize FD/RoaQ)
summary(full22.mod <- aov(full~trt, subdatFD22)) #summary model
full22.tuk <- generateTukeyLabel(full22.mod, subdatFD22$full) #tukey and labels
FD.wy.22 <- subdatFD22 %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq22.wy$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=subdatFD22, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels






#### 2023 ####
summary(leafn23.mod <- aov(leafn~trt, subdat23))
leafn23.tuk <- generateTukeyLabel(leafn23.mod, subdat23$leafn)
leafn.wy.23 <- CWM_trait_plot(subdat23,subdat23$leafn,"Leaf N") +
  geom_point(aes(y=quantile(leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
  #  geom_point(aes(y=quantile(cwm23.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl23.mod <- aov(srl~trt, subdat23))
srl23.tuk <- generateTukeyLabel(srl23.mod, subdat23$srl)
srl.wy.23 <- CWM_trait_plot(subdat23,subdat23$srl,"Specific root length") +
  geom_point(aes(y=quantile(srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = srl23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

# veg.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$veg,"Vegetative spread potential") +
#   ylim(0,1)
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD23$veg <- normalize(subdatFD23$veg) #(normalize FD/RoaQ)
summary(veg23.mod <- aov(srl~trt, subdatFD23)) #summary model
veg23.tuk <- generateTukeyLabel(veg23.mod, subdatFD23$veg) #tukey and labels
veg.wy.23 <- CWM_trait_FD_plot(subdatFD23,subdatFD23$veg,"Vegetative spread potential") +
  geom_point(aes(y=max(veg),x="ir"),data=subdatFD23, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = veg23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=2, col="black") #add tukey labels

summary(ldmc23.mod <- aov(ldmc~trt, subdat23))
ldmc23.tuk <- generateTukeyLabel(ldmc23.mod, subdat23$ldmc)
ldmc.wy.23 <- CWM_trait_plot(subdat23,subdat23$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = ldmc23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lop23.mod <- aov(lop~trt, subdat23))
lop23.tuk <- generateTukeyLabel(lop23.mod, subdat23$lop)
lop.wy.23 <- CWM_trait_plot(subdat23,subdat23$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = lop23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2, vjust=-1.75, col="black") #add tukey labels

# rootdiam.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD23$rootdiam <- normalize(subdatFD23$rootdiam) #(normalize FD/RoaQ)
summary(rd23.mod <- aov(rootdiam~trt, subdatFD23)) #summary model
rd23.tuk <- generateTukeyLabel(rd23.mod, subdatFD23$rootdiam) #tukey and labels
rootdiam.wy.23 <- CWM_trait_FD_plot(subdatFD23,subdatFD23$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=subdatFD23, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

subdatFD23$full <- normalize(subdatFD23$full) #(normalize FD/RoaQ)
summary(full23.mod <- aov(full~trt, subdatFD23)) #summary model
full23.tuk <- generateTukeyLabel(full23.mod, subdatFD23$full) #tukey and labels
FD.wy.23 <- subdatFD23 %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq23.wy$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=subdatFD23, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels


## Produce a .tiff file of our plots together! (if I want these pots, rename things below)
library(patchwork)
# export 21
tiff("figures/cwm WY/traits_2021(plot).tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.wy.21 + srl.wy.21 + veg.wy.21) / (ldmc.wy.21 + lop.wy.21 + rootdiam.wy.21)
dev.off()
tiff("figures/cwm WY/FD_2021.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.wy.21
dev.off()
# export 22
tiff("figures/cwm WY/traits_2022(plot).tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.wy.22 + srl.wy.22 + veg.wy.22) / (ldmc.wy.22 + lop.wy.22 + rootdiam.wy.22)
dev.off()
tiff("figures/cwm WY/FD_2022.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.wy.22
dev.off()
# export 23
tiff("figures/cwm WY/traits_2023(plot).tiff", res=400, height = 6,width =8.5, "in",compression = "lzw")
(leafn.wy.23 + srl.wy.23 + veg.wy.23) / (ldmc.wy.23 + lop.wy.23 + rootdiam.wy.23)
dev.off()
tiff("figures/cwm WY/FD_2023.tiff", res=400, height = 4,width =4, "in",compression = "lzw")
FD.wy.23
dev.off()

## figure for short report
#21
tiff("figures/cwm wy/alltargets_2021(plot).tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.wy.21 + srl.wy.21 + veg.wy.21) / (ldmc.wy.21 + lop.wy.21 + rootdiam.wy.21)) | (FD.wy.21)) +
  plot_layout(widths = c(2,1)) +
  plot_annotation(title = 'WY 2021')
dev.off()
#22
tiff("figures/cwm wy/alltargets_2022(plot).tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.wy.22 + srl.wy.22 + veg.wy.22) / (ldmc.wy.22 + lop.wy.22 + rootdiam.wy.22)) | (FD.wy.22)) +
  plot_layout(widths = c(2,1)) +
  plot_annotation(title = 'WY 2022')
dev.off()
#23
tiff("figures/cwm wy/alltargets_2023(plot).tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.wy.23 + srl.wy.23 + veg.wy.23) / (ldmc.wy.23 + lop.wy.23 + rootdiam.wy.23)) | (FD.wy.23)) +
  plot_layout(widths = c(2,1)) +
  plot_annotation(title = 'WY 2023')
dev.off()