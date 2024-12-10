#### First, for supplemental, create summary PCA of all species by 
#### their traits and PCA of communities by seeding treatment 
#### Next, make NMDS of communities by treatment for conceptual 
#### diagram.
#### All figures use predicted/ model-derived species abundances 

#library loading
library(tidyverse)
library(vegan)
library(RColorBrewer)

#### WY
### Create PCA of species traits
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
subtraits <- traits.wy %>% select(leafn,srl,rootdiam,lop,ldmc,veg)
#Subset monocots and dicots
grams <- subset(traits.wy, graminoid==1)
forbs <- subset(traits.wy, graminoid==0)

## Make pca
pca.wy <- princomp(subtraits, cor=TRUE)
summary(pca.wy)
biplot(pca.wy)
traits.wy$pc1 <- pca.wy$scores[,1]
traits.wy$pc2 <- pca.wy$scores[,2]

library("factoextra")
wy.pca.plot <- fviz_pca_biplot(pca.wy, geom = c("point","text"), ggtheme = theme_minimal(), 
                               alpha.ind=0.8, col.var = "black", repel=TRUE, labelsize = 3,
                               habillage =traits.wy$graminoid, title = "Wyoming")+
  guides(shape="none", fill="none")+
  scale_color_manual(values=c("tan4","#009E73"), labels=c("Forb","Grass"))+
  labs(col="Lifeform")

#### CA
### Create PCA of species traits
# CSV of species-trait combinations (for OG 32)
traits.ca <- read.csv("data/trait_data/annualgrass.csv", header=TRUE, row.names=1)
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
rownames(traits.ca)<-traits.ca$Code

### pca
pca.ca <- princomp(traits.ca[,c(2:4,7,9,10)], cor=TRUE)
summary(pca.ca)
biplot(pca.ca)
traits.ca$pc1 <- pca$scores[,1]
traits.ca$pc2 <- pca$scores[,2]

ca.pca.plot <- fviz_pca_biplot(pca.ca, geom = c("point","text"), ggtheme = theme_minimal(), 
                               alpha.ind=0.8, col.var = "black", repel=TRUE, labelsize = 3,
                               habillage =traits.ca$graminoid, title = "California")+
  guides(shape="none", fill="none")+
  scale_color_manual(values=c("tan4","#009E73"), labels=c("Forb","Grass"))+
  scale_shape_manual(values=c("circle","triangle"), labels=c("Forb","Grass"))+
  labs(col="Lifeform")

### combine species PCA's
library(ggpubr)
#export for report
tiff("figures/sp_pca_plots.tiff", res=400, height = 6,width =6, "in",compression = "lzw")
ggarrange(wy.pca.plot,ca.pca.plot, nrow=2,common.legend = T)
dev.off()


#### WY
### continue wrangling trait data to use to create CWM pca
#set color scheme
colors = brewer.pal(12,"Paired")
traits.wy$cols = c(colorRampPalette(colors)(nrow(traits.wy)))
# Arrange trait matrix alphabetically
trait.matrix.wy <- traits.wy[order(rownames(traits.wy)),]
#Select colons of interest
trait.matrix.wy <- as.matrix(trait.matrix.wy[,2:11])

### make model-derived CWM PCA
## pretreatment/ seeding probability 
preds.wy <- read.csv("data/allplot.assemblages.csv") #data
preds.wy <- preds.wy %>% arrange(trt,block)
comms_p.wy <- labdsv::matrify(data.frame(preds.wy$trt,preds.wy$species,preds.wy$prob))
comms_p.wy <- comms_p.wy[,order(colnames(comms_p.wy))]
og25 <- unique(colnames(comms_p.wy)) #make list of 25 native species
## Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
cwm_p.wy <- FD::functcomp(as.matrix(trait.matrix.wy), as.matrix(comms_p.wy), bin.num=c("graminoid","veg","c4"))
#Define block by extracting the numeric from the cwm rownames
cwm_p.wy$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(cwm_p.wy)))
#Define trt by extracting the subplot from the cwm rownames
cwm_p.wy$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(cwm_p.wy))))
# Define drought treatment at block level
covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
cwm_p.wy <- cwm_p.wy %>% mutate(drought = case_when(block %in% covered ~ "drt",
                                                    !block %in% covered ~ "cntl")) 
cwm_p.wy 
## subset for just CWM per trt
subcwm_p.wy <- cwm_p.wy %>% select(-sla,-rdmc,-c4,-drought,-block,-trt) #[-c(7,10:13),]
comms_p.wy <- comms_p.wy #community/species data

## create PCA by year for each site
# pcawy <- princomp(cwm_p.wy)
# biplot(pcawy)
pca_result <- prcomp(subcwm_p.wy) #PCA traits
pca_result2 <- prcomp(comms_p.wy) #PCA species

# Extract scores for the PCA plot
pca_scores <- as.data.frame(pca_result$x)
pca_scores2 <- as.data.frame(pca_result2$x)

#Define block by extracting the numeric from the cwm rownames
pca_scores$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(pca_scores)))
pca_scores2$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(pca_scores2)))
#Define trt by extracting the subplot from the cwm rownames
pca_scores$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(pca_scores))))
pca_scores2$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(pca_scores2))))
# Define drought treatment at block level
# covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
# pca_result <- pca_result %>% mutate(drought = case_when(block %in% covered ~ "drt",
 
###don't make these plots ggplot does not make ellipses correctly
# # Step 3: Plot PCA with ggplot2, adding ellipses by treatment group
# ggplot(pca_scores, aes(x = PC1, y = PC2, color = trt)) +
#   geom_point(alpha = 0.6) +  # Points for each plot
#   stat_ellipse(aes(fill = trt), geom = "polygon", alpha = 0.2) +  # Ellipses for each treatment
#   labs(x = "PC1", y = "PC2", title = "PCA of Functional Traits by Seeding Treatment") +
#   theme_minimal() +
#   scale_color_viridis_d() +
#   scale_fill_viridis_d() +
#   theme(legend.position = "right")
# 
# ggplot(pca_scores2, aes(x = PC1, y = PC2, color = trt)) +
#   geom_point(alpha = 0.6) +  # Points for each plot
#   stat_ellipse(aes(fill = trt), geom = "polygon", alpha = 0.2) +  # Ellipses for each treatment
#   labs(x = "PC1", y = "PC2", title = "PCA of Species by Seeding Treatment") +
#   theme_minimal() +
#   scale_color_viridis_d() +
#   scale_fill_viridis_d() +
#   theme(legend.position = "right")

### NMDS of communities;
#run Bray-Curtis with vegan::vegdist() 
testbray.wy <- vegdist(as.matrix(comms_p.wy), method = "bray")
## make view-able matrix + subset so 2020 comms are rows and 2024 comms are columns 
testbray.mat.wy <- as.matrix(testbray.wy)
nmdstest.wy <- metaMDS(testbray.mat.wy, distance = "bray")
# plot(nmdstest)
# orditorp(nmdstest, "sites")
ordiplot(testbray.wy, type = "points", main = "ellipses")
points(testbray.wy, display = "sites", col = as.numeric(groupswy), pch = 16)
ordiellipse(testbray.wy, groupswy, kind = "sd", conf = 0.95, display = "sites")
groupswy<- as.factor(nmds_scores.wy$trt)
ordiellipse(nmdstest.wy, groups = trt, draw = "polygon", lty = 1, col = "grey90")
nmds_scores.wy <- as.data.frame(scores(nmdstest.wy))


colnames(ellipses) <- c("NMDS1", "NMDS2", "trt")

# Plot
ggplot(nmds_scores.wy, aes(x = NMDS1, y = NMDS2, color = trt)) +
  geom_point(alpha = 0.6) +
  geom_path(data = ellipses, aes(x = NMDS1, y = NMDS2, color = trt), size = 1) +
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS with Ellipses") +
  theme_minimal()
##whatever ^


#Define block by extracting the numeric from the cwm rownames
nmds_scores.wy$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(nmds_scores.wy)))
#Define trt by extracting the subplot from the cwm rownames
nmds_scores.wy$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(nmds_scores.wy))))

# Create plot with ggplot2
ggplot(nmds_scores.wy, aes(x = NMDS1, y = NMDS2, color = trt)) +
  geom_point(alpha = 0.6) +  # Points for each plot
  stat_ellipse(aes(fill = trt), type = "norm", level = 0.95, alpha = 0.2) +  # Ellipses for each treatment
  labs(x = "NMDS1", y = "NMDS2", title = "PCA of Species by Seeding Treatment") +
  theme_minimal() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.position = "right")



#### CA
### continue wrangling trait data to use to create CWM pca
#set color scheme
colors = brewer.pal(12,"Paired")
traits.ca$cols = c(colorRampPalette(colors)(nrow(traits.ca)))
## Subset monocots and dicots
grams <- subset(traits.ca, graminoid==1)
forbs <- subset(traits.ca, graminoid==0)
traits.ca$graminoid<-as.factor(traits.ca$graminoid)
# Arrange trait matrix alphabetically
trait.matrix.ca <- traits.ca[order(rownames(traits.ca)),]
#Select columns of interest, make into matrix
test <- data.frame(trait.matrix.ca[,2:8], row.names = trait.matrix.ca[,1])
rownames(test)
trait.matrix.ca <- as.matrix(test)
trait.matrix.ca <- trait.matrix.ca[order(rownames(trait.matrix.ca)),]

### make model-derived CWM PCA
## pretreatment/ seeding probability 
preds.ca <- read.csv("data/calgrass.allplot.assemblages.csv") #data
preds.ca$trt.b <- paste(preds.ca$trt, preds.ca$block)
preds.ca <- preds.ca %>% arrange(trt.b,block)
comms_p.ca <- labdsv::matrify(data.frame(preds.ca$trt.b,preds.ca$species,preds.ca$prob))
comms_p.ca <- comms_p.ca[,order(colnames(comms_p.ca))]
#comms_p.ca <- comms_p.ca %>% select(-CAME)

# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
trait.matrix.ca.pred <- traits.ca[order(rownames(traits.ca)),]
test <- data.frame(trait.matrix.ca.pred[,c(2:4,7,9,10:12)], row.names = trait.matrix.ca.pred[,1])
rownames(test)
trait.matrix.ca.pred <- as.matrix(test)
trait.matrix.ca.pred <- trait.matrix.ca.pred[order(rownames(trait.matrix.ca.pred)),]
cwm_p.ca <- FD::functcomp(as.matrix(trait.matrix.ca.pred), as.matrix(comms_p.ca), bin.num=c("graminoid"))
#Define block by extracting the numeric from the cwm rownames
cwm_p.ca$block <- as.factor(sub("^[a-z]+ (\\d+)$", "\\1", rownames(cwm_p.ca)))
#Define seeding trt by extracting the letters from the cwm rownames
cwm_p.ca$trt <- as.factor(sub("(^[a-z]+) \\d+$", "\\1", rownames(cwm_p.ca)))
cwm_p.ca <- cwm_p.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand 
# # Define drought treatment at block level
# block.water <- comp.ca %>% select(c(block,trt,water))# get data from OG dataframe
# block.water$trt <- tolower(block.water$trt) #make lowercase to match
# block.water$trt<-replace(block.water$trt,block.water$trt=="r","rand") #make r and rand match
# block.water <- block.water[c(1:210),] #repeating 3 times(?) 
# cwm_p.ca <- merge(cwm_p.ca, block.water, by = c("block","trt"), all.x=T) #merge

##from cwm_trait_ca
subcwm_p.ca <- cwm_p.ca %>% select(-RTD,-trt,-block)
comms_p.ca <- comms_p.ca

# ## create PCA by year for each site
# pcaca <- princomp(ogcadat)
# biplot(pcaca)
pca_result <- prcomp(subcwm_p.ca) #traits PCA
pca_result2 <- prcomp(comms_p.ca) #species PCA

# Extract scores for the PCA plot
pca_scores <- as.data.frame(pca_result$x)
pca_scores2 <- as.data.frame(pca_result2$x)

#Define block by extracting the numeric from the cwm rownames
pca_scores$block <- as.factor(sub("^[a-z]+ (\\d+)$", "\\1", rownames(pca_scores)))
pca_scores2$block <- as.factor(sub("^[a-z]+ (\\d+)$", "\\1", rownames(pca_scores2)))
#Define seeding trt by extracting the letters from the cwm rownames
pca_scores$trt <- as.factor(sub("(^[a-z]+) \\d+$", "\\1", rownames(pca_scores)))
pca_scores <- pca_scores %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand 
pca_scores2$trt <- as.factor(sub("(^[a-z]+) \\d+$", "\\1", rownames(pca_scores2)))
pca_scores2 <- pca_scores2 %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand 
# covered <- as.character(c(3,4,6,7,9,11,14,15,18,20,22,24,26,29,30,32,34,35,36,39,41,45,47,50,53,54,56,58,59,61,62,63))
# pca_result <- pca_result %>% mutate(drought = case_when(block %in% covered ~ "drt",

# # not using now.
# # Step 3: Plot PCA with ggplot2, adding ellipses by treatment group
# ggplot(pca_scores, aes(x = PC1, y = PC2, color = trt)) +
#   geom_point(alpha = 0.6) +  # Points for each plot
#   stat_ellipse(aes(fill = trt), geom = "polygon", alpha = 0.2) +  # Ellipses for each treatment
#   labs(x = "PC1", y = "PC2", title = "PCA of Functional Traits by Seeding Treatment") +
#   theme_minimal() +
#   scale_color_viridis_d() +
#   scale_fill_viridis_d() +
#   theme(legend.position = "right")
# 
# ggplot(pca_scores2, aes(x = PC1, y = PC2, color = trt)) +
#   geom_point(alpha = 0.6) +  # Points for each plot
#   stat_ellipse(aes(fill = trt), geom = "polygon", alpha = 0.2) +  # Ellipses for each treatment
#   labs(x = "PC1", y = "PC2", title = "PCA of Species by Seeding Treatment") +
#   theme_minimal() +
#   scale_color_viridis_d() +
#   scale_fill_viridis_d() +
#   theme(legend.position = "right")


######## NMDS 
#run Bray-Curtis with vegan::vegdist() 
testbray.ca <- vegdist(as.matrix(comms_p.ca), method = "bray")
testbray.wy <- vegdist(as.matrix(comms_p.wy), method = "bray")

## make view-able matrix + subset so 2020 comms are rows and 2024 comms are columns 
testbray.mat.ca <- as.matrix(testbray.ca)
testbray.mat.wy <- as.matrix(testbray.wy)

nmdstest.ca <- metaMDS(testbray.mat.ca, distance = "bray")
nmdstest.wy <- metaMDS(testbray.mat.wy, distance = "bray")

nmds_scores.ca <- as.data.frame(scores(nmdstest.ca))
nmds_scores.wy <- as.data.frame(scores(nmdstest.wy))

#Define block by extracting the numeric from the cwm rownames
nmds_scores.ca$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(nmds_scores.ca)))
nmds_scores.wy$block <- as.factor(sub(".*\\s(\\d+)$", "\\1", rownames(nmds_scores.wy)))
#Define trt by extracting the subplot from the cwm rownames
nmds_scores.ca$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(nmds_scores.ca))))
nmds_scores.wy$trt <- as.factor(as.character(sub("(.*\\s)\\d+$", "\\1", rownames(nmds_scores.wy))))

# Create plot with ggplot2
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = as.factor(trt))) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(type = "t", level = 0.95) +  # Add ellipses for confidence interval
  scale_color_viridis_d(name = "Group") +
  theme_classic() +
  labs(
    title = "NMDS of Plant Communities",
    x = "NMDS1",
    y = "NMDS2"
  )

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = trt)) +
  geom_point(alpha = 0.6) +  # Points for each plot
  stat_ellipse(aes(fill = trt), geom = "polygon", alpha = 0.2) +  # Ellipses for each treatment
  labs(x = "NMDS1", y = "NMDS2", title = "PCA of Species by Seeding Treatment") +
  theme_minimal() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.position = "right")

ggplot(nmds_scores.wy, aes(x = NMDS1, y = NMDS2, color = trt)) +
  geom_point(alpha = 0.6) +  # Points for each plot
  stat_ellipse(aes(fill = trt), geom = "polygon", alpha = 0.2) +  # Ellipses for each treatment
  labs(x = "NMDS1", y = "NMDS2", title = "PCA of Species by Seeding Treatment") +
  theme_minimal() +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.position = "right")

### create plots
library(ggordiplots)
ordiplot_output<-gg_ordiplot(nmdstest.wy, nmds_scores.wy$trt, choices=c(1,2), ellipse = T)#+
  scale_color_viridis_d()
gg_ordiplot(nmdstest.ca, trtca)
trtca <- nmds_scores.ca$trt

df_ellipse <- ordiplot_output$df_ellipse
df_ellipse$Group <- as.factor(df_ellipse$Group)
dt_ellipse <- df_ellipse %>% filter(Group == "dt ")
fd_ellipse <- subset(df_ellipse, Group == "fd ", droplevels= TRUE)
ir_ellipse <- subset(df_ellipse, Group == "ir ", droplevels= TRUE)
r_ellipse <- subset(df_ellipse, Group == "r ", droplevels= TRUE)

df_pt <- ordiplot_output$df_ord
df_pt$Group <- as.factor(df_pt$Group)
dt_pt <- df_pt %>% filter(Group == "dt ")
fd_pt <- subset(df_pt, Group == "fd ", droplevels= TRUE)
ir_pt <- subset(df_pt, Group == "ir ", droplevels= TRUE)
r_pt <- subset(df_pt, Group == "r ", droplevels= TRUE)

nmdsplot.wy <- ggplot(nmds_scores.wy, aes(x = NMDS1, y = NMDS2, color = trt)) +
  geom_point(data=dt_pt, aes(x=x, y=y), color="#3E4A89", alpha = 0.6)+
  geom_point(data=fd_pt, aes(x=x, y=y), color="#35B779", alpha = 0.6)+
  geom_point(data=ir_pt, aes(x=x, y=y), color="#FDE725", alpha = 0.6)+
  geom_point(data=r_pt, aes(x=x, y=y), color="#440154", alpha = 0.6)+
  geom_path(data=dt_ellipse, aes(x=x, y=y), color="#3E4A89", linetype="solid", size=1) +
  geom_path(data=fd_ellipse, aes(x=x, y=y), color="#35B779", linetype="solid", size=1) +
  geom_path(data=ir_ellipse, aes(x=x, y=y), color="#FDE725", linetype="solid", size=1) +
  geom_path(data=r_ellipse, aes(x=x, y=y), color="#440154", linetype="solid", size=1)+
  #geom_polygon(data = dt_ellipse, aes(x=x, y=y), fill="#3E4A89",group=trt, alpha = 0.3)+
  # geom_polygon(data = fd_ellipse, aes(x=x, y=y), fill="#35B779", alpha = 0.2)+
  # geom_polygon(data = ir_ellipse, aes(x=x, y=y), fill="#FDE725", alpha = 0.6)+
  # geom_polygon(data = r_ellipse, aes(x=x, y=y), fill="#440154", alpha = 0.2)+
  # scale_fill_viridis_d()+
  #geom_polygon(data = dt_ellipse, aes(x = x, y = y, fill = Group, group = Group), alpha = 0.2) +
  labs(x = "NMDS1", y = "NMDS2", title = "PCA of communities by seeding treatment") +
  theme_minimal()


ordiplot_output.ca<-gg_ordiplot(nmdstest.ca, nmds_scores.ca$trt, choices=c(1,2), ellipse = T)#+

df_ellipse.ca <- ordiplot_output.ca$df_ellipse
df_ellipse.ca$Group <- as.factor(df_ellipse.ca$Group)
dt_ellipse.ca <- df_ellipse.ca %>% filter(Group == "dt ")
fd_ellipse.ca <- subset(df_ellipse.ca, Group == "fd ", droplevels= TRUE)
ir_ellipse.ca <- subset(df_ellipse.ca, Group == "ir ", droplevels= TRUE)
r_ellipse.ca <- subset(df_ellipse.ca, Group == "r ", droplevels= TRUE)

df_pt.ca <- ordiplot_output.ca$df_ord
df_pt.ca$Group <- as.factor(df_pt.ca$Group)
dt_pt.ca <- df_pt.ca %>% filter(Group == "dt ")
fd_pt.ca <- subset(df_pt.ca, Group == "fd ", droplevels= TRUE)
ir_pt.ca <- subset(df_pt.ca, Group == "ir ", droplevels= TRUE)
r_pt.ca <- subset(df_pt.ca, Group == "r ", droplevels= TRUE)

tester <- ordiplot_output.ca$df_mean.ord
tester$Group <- as.factor(tester$Group)
tester2 <- tester %>% filter(Group == "dt ")

nmdsplot.ca <- ggplot(nmds_scores.ca, aes(x = NMDS1, y = NMDS2, color = trt)) +
  geom_point(data=dt_pt.ca, aes(x=x, y=y), color="#3E4A89", alpha = 0.6)+
  geom_point(data=fd_pt.ca, aes(x=x, y=y), color="#35B779", alpha = 0.6)+
  geom_point(data=ir_pt.ca, aes(x=x, y=y), color="#FDE725", alpha = 0.6)+
  geom_point(data=r_pt.ca, aes(x=x, y=y), color="#440154", alpha = 0.6)+
  geom_path(data=dt_ellipse.ca, aes(x=x, y=y), color="#3E4A89", linetype="solid", size=1) +
  geom_path(data=fd_ellipse.ca, aes(x=x, y=y), color="#35B779", linetype="solid", size=1) +
  geom_path(data=ir_ellipse.ca, aes(x=x, y=y), color="#FDE725", linetype="solid", size=1) +
  geom_path(data=r_ellipse.ca, aes(x=x, y=y), color="#440154", linetype="solid", size=1)+
  #geom_circle(data = tester2, aes(x=x, y=y),alpha = 0.3)+
  #geom_polygon(data = dt_ellipse.ca, aes(x=x, y=y, fill=Group),alpha = 0.3)+
  # geom_polygon(data = fd_ellipse, aes(x=x, y=y), fill="#35B779", alpha = 0.2)+
  # geom_polygon(data = ir_ellipse, aes(x=x, y=y), fill="#FDE725", alpha = 0.6)+
  # geom_polygon(data = r_ellipse, aes(x=x, y=y), fill="#440154", alpha = 0.2)+
  # scale_fill_viridis_d()+
  #geom_polygon(data = dt_ellipse, aes(x = x, y = y, fill = Group, group = Group), alpha = 0.2) +
  labs(x = "NMDS1", y = "NMDS2", title = "PCA of communities by seeding treatment") +
  theme_minimal()

library(ggpubr)
ggarrange(nmdsplot.wy,nmdsplot.ca)
tiff("figures/NMDSplots.tiff", res=400, height = 2.5,width =5, "in",compression = "lzw")
ggarrange(nmdsplot.wy,nmdsplot.ca, ncol=2, common.legend = T)
dev.off()