#### General assessment of CWM in different seeding treatment in WY
#### To see how CWM's differed by seeding trt, the following code 
#### creates violin plots and test for differences using ANOVA 
#### (as well as +post-hoc Tukey and plotting) annually.
#### These analyses and figures were used for visualization and initial
#### assessment, but none of these analyses or figures are included in
#### the final manuscript. 

#library loading
library(tidyverse)

#### Create functions
## Functions needed to normalizing RoaQ, making post-hoc letters, 
## and creating figures
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
dat <- read.csv("data/cwm_wy.csv") #cwm data
dat$trt <- factor(dat$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
datFD <- read.csv("data/cwm_raoq_wy.csv") #cwm RoaQ data
datFD$trt <- factor(datFD$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
#trait data
traits.wy <- read.csv("data/mixedgrass.csv", header=TRUE, row.names=1) # CSV of species-trait combinations (for OG 25)
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

## figure for short report
#21
tiff("figures/cwm wy/alltargets_2021.tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.wy.21 + srl.wy.21 + veg.wy.21) / (ldmc.wy.21 + lop.wy.21 + rootdiam.wy.21)) | (FD.wy.21)) +
  plot_layout(widths = c(2,1)) + 
  plot_annotation(title = 'WY 2021')
dev.off()
#22
tiff("figures/cwm wy/alltargets_2022.tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.wy.22 + srl.wy.22 + veg.wy.22) / (ldmc.wy.22 + lop.wy.22 + rootdiam.wy.22)) | (FD.wy.22)) +
  plot_layout(widths = c(2,1)) + 
  plot_annotation(title = 'WY 2022')
dev.off()
#23
tiff("figures/cwm wy/alltargets_2023.tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.wy.23 + srl.wy.23 + veg.wy.23) / (ldmc.wy.23 + lop.wy.23 + rootdiam.wy.23)) | (FD.wy.23)) +
  plot_layout(widths = c(2,1)) + 
  plot_annotation(title = 'WY 2023')
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


## All years combined
summary(leafn.mod <- aov(leafn~trt, dat))
leafn.tuk <- generateTukeyLabel(leafn.mod, dat$leafn)
leafn.wy <- CWM_trait_plot(dat,dat$leafn,"Leaf N") +
  geom_point(aes(y=quantile(leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
  #  geom_point(aes(y=quantile(cwm21.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl.mod <- aov(srl~trt, dat))
srl.tuk <- generateTukeyLabel(srl.mod, dat$srl)
srl.wy <- CWM_trait_plot(dat,dat$srl,"Specific root length") +
  geom_point(aes(y=quantile(srl,.7557),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = srl.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

# veg.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$veg,"Vegetative spread potential") +
#   ylim(0,1)
#plot veg as Raoq instead (normalize FD/RoaQ)
datFD$veg <- normalize(datFD$veg) #(normalize FD/RoaQ)
summary(veg.mod <- aov(veg~trt, datFD)) #summary model
veg.tuk <- generateTukeyLabel(veg.mod, datFD$veg) #tukey and labels
veg.wy <- CWM_trait_FD_plot(datFD,datFD$veg,"Vegetative spread potential") +
  geom_point(aes(y=max(veg),x="ir"),data=datFD, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = veg.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=2, col="black") #add tukey labels

summary(ldmc.mod <- aov(ldmc~trt, dat))
ldmc.tuk <- generateTukeyLabel(ldmc.mod, dat$ldmc)
ldmc.wy <- CWM_trait_plot(dat,dat$ldmc,"Leaf dry matter content") +
  geom_point(aes(y=quantile(ldmc,.75),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = ldmc.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lop.mod <- aov(lop~trt, dat))
lop.tuk <- generateTukeyLabel(lop.mod, dat$lop)
lop.wy <- CWM_trait_plot(dat,dat$lop,"Leaf osmotic potential") +
  geom_point(aes(y=quantile(lop,.25),x="dt"),data=traits.wy, col= "red", shape=18, size=3.5) +
  geom_text(data = lop.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2, vjust=-1.75, col="black") #add tukey labels

# rootdiam.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
datFD$rootdiam <- normalize(datFD$rootdiam) #(normalize FD/RoaQ)
summary(rd.mod <- aov(rootdiam~trt, datFD)) #summary model
rd.tuk <- generateTukeyLabel(rd.mod, datFD$rootdiam) #tukey and labels
rootdiam.wy <- CWM_trait_FD_plot(datFD,datFD$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=datFD, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

datFD$full <- normalize(datFD$full) #(normalize FD/RoaQ)
summary(full.mod <- aov(full~trt, datFD)) #summary model
full.tuk <- generateTukeyLabel(full.mod, datFD$full) #tukey and labels
FD.wy <- datFD %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq21.wy$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=datFD, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels


#figure together
png("figures/cwm wy/alltargets.png", res=400, height = 6,width =12, "in")#,compression = "lzw")
(((leafn.wy + srl.wy + veg.wy) / (ldmc.wy + lop.wy + rootdiam.wy)) | (legend)/(FD.wy)) +
  plot_layout(widths = c(2.25,.75)) + 
  plot_annotation(title = 'WY (all years)')
dev.off()