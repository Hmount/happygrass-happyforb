#### CWM ANOVA's and figures for CA site
#### How to CWM's differ by trt?
#### How do CWM's maintain/change through time? (not on here right now)

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
dat <- read.csv("data/cwm_ca.csv") #cwm data
dat$trt <- factor(dat$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
datFD <- read.csv("data/cwm_raoq_ca.csv") #cwm RoaQ data
datFD$trt <- factor(datFD$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
#trait data
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

## subset data by year for figures
subdat21 <- dat %>% filter(year=="2021")
subdatFD21 <- datFD %>% filter(year=="2021")
subdat22 <- dat %>% filter(year=="2022")
subdatFD22 <- datFD %>% filter(year=="2022")
subdat23 <- dat %>% filter(year=="2023")
subdatFD23 <- datFD %>% filter(year=="2023")

#### Did traits of seeding treatments differ? (ANOVA, post-hoc Tukey, and plotting)
#### 2021 ####
summary(leafn21.mod <- aov(N~trt, subdat21))
leafn21.tuk <- generateTukeyLabel(leafn21.mod, subdat21$N)
leafn.ca.21 <- CWM_trait_plot(subdat21,subdat21$N,"Leaf N") +
  geom_point(aes(y=quantile(N,.33),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
#  geom_point(aes(y=quantile(cwm21.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl21.mod <- aov(SRL~trt, subdat21))
srl21.tuk <- generateTukeyLabel(srl21.mod, subdat21$SRL)
srl.ca.21 <- CWM_trait_plot(subdat21,subdat21$SRL,"Specific root length") +
  geom_point(aes(y=quantile(SRL,.67),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = srl21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=-1.75, col="black") #add tukey labels

summary(rmf21.mod <- aov(RMF~trt, subdat21))
rmf21.tuk <- generateTukeyLabel(rmf21.mod, subdat21$RMF)
rmf.ca.21 <- CWM_trait_plot(subdat21,subdat21$RMF,"Root mass fraction") +
  geom_point(aes(y=quantile(RMF,.33),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = rmf21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lma21.mod <- aov(LMA~trt, subdat21))
lma21.tuk <- generateTukeyLabel(lma21.mod, subdat21$LMA)
lma.ca.21 <- CWM_trait_plot(subdat21,subdat21$LMA,"Leaf mass per area") +
  geom_point(aes(y=quantile(LMA,.67),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = lma21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(seedmass21.mod <- aov(seed.mass~trt, subdat21))
seedmass21.tuk <- generateTukeyLabel(seedmass21.mod, subdat21$seed.mass)
seedmass.ca.21 <- CWM_trait_plot(subdat21,subdat21$seed.mass,"Seed mass") +
  geom_point(aes(y=quantile(seed.mass,.67),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = seedmass21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

# rootdiam.wy.21 <- CWM_trait_plot(cwm21.wy,cwm21.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD21$rootdiam <- normalize(subdatFD21$rootdiam) #(normalize FD/RoaQ)
summary(rd21.mod <- aov(rootdiam~trt, subdatFD21)) #summary model
rd21.tuk <- generateTukeyLabel(rd21.mod, subdatFD21$rootdiam) #tukey and labels
rootdiam.ca.21 <- CWM_trait_FD_plot(subdatFD21,subdatFD21$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=subdatFD21, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

subdatFD21$full <- normalize(subdatFD21$full) #(normalize FD/RoaQ)
summary(full21.mod <- aov(full~trt, subdatFD21)) #summary model
full21.tuk <- generateTukeyLabel(full21.mod, subdatFD21$full) #tukey and labels
FD.ca.21 <- subdatFD21 %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq21.ca$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=subdatFD21, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full21.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels



#### 2022 ####
summary(leafn22.mod <- aov(N~trt, subdat22))
leafn22.tuk <- generateTukeyLabel(leafn22.mod, subdat22$N)
leafn.ca.22 <- CWM_trait_plot(subdat22,subdat22$N,"Leaf N") +
  geom_point(aes(y=quantile(N,.33),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
  #  geom_point(aes(y=quantile(cwm22.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl22.mod <- aov(SRL~trt, subdat22))
srl22.tuk <- generateTukeyLabel(srl22.mod, subdat22$SRL)
srl.ca.22 <- CWM_trait_plot(subdat22,subdat22$SRL,"Specific root length") +
  geom_point(aes(y=quantile(SRL,.67),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = srl22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=-1.75, col="black") #add tukey labels

summary(rmf22.mod <- aov(RMF~trt, subdat22))
rmf22.tuk <- generateTukeyLabel(rmf22.mod, subdat22$RMF)
rmf.ca.22 <- CWM_trait_plot(subdat22,subdat22$RMF,"Root mass fraction") +
  geom_point(aes(y=quantile(RMF,.33),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = rmf22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lma22.mod <- aov(LMA~trt, subdat22))
lma22.tuk <- generateTukeyLabel(lma22.mod, subdat22$LMA)
lma.ca.22 <- CWM_trait_plot(subdat22,subdat22$LMA,"Leaf mass per area") +
  geom_point(aes(y=quantile(LMA,.67),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = lma22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(seedmass22.mod <- aov(seed.mass~trt, subdat22))
seedmass22.tuk <- generateTukeyLabel(seedmass22.mod, subdat22$seed.mass)
seedmass.ca.22 <- CWM_trait_plot(subdat22,subdat22$seed.mass,"Seed mass") +
  geom_point(aes(y=quantile(seed.mass,.67),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = seedmass22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

# rootdiam.wy.22 <- CWM_trait_plot(cwm22.wy,cwm22.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD22$rootdiam <- normalize(subdatFD22$rootdiam) #(normalize FD/RoaQ)
summary(rd22.mod <- aov(rootdiam~trt, subdatFD22)) #summary model
rd22.tuk <- generateTukeyLabel(rd22.mod, subdatFD22$rootdiam) #tukey and labels
rootdiam.ca.22 <- CWM_trait_FD_plot(subdatFD22,subdatFD22$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=subdatFD22, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

subdatFD22$full <- normalize(subdatFD22$full) #(normalize FD/RoaQ)
summary(full22.mod <- aov(full~trt, subdatFD22)) #summary model
full22.tuk <- generateTukeyLabel(full22.mod, subdatFD22$full) #tukey and labels
FD.ca.22 <- subdatFD22 %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq22.ca$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=subdatFD22, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full22.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels



#### 2023 ####
summary(leafn23.mod <- aov(N~trt, subdat23))
leafn23.tuk <- generateTukeyLabel(leafn23.mod, subdat23$N)
leafn.ca.23 <- CWM_trait_plot(subdat23,subdat23$N,"Leaf N") +
  geom_point(aes(y=quantile(N,.33),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) + # Specifying the target object (red dot).
  #  geom_point(aes(y=quantile(cwm23.wy$leafn,.25),x="ir"), data=traits.wy, col= "red", shape=18, size=3.5) # target could be specified per year, but probably makes less sense.
  geom_text(data = leafn23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(srl23.mod <- aov(SRL~trt, subdat23))
srl23.tuk <- generateTukeyLabel(srl23.mod, subdat23$SRL)
srl.ca.23 <- CWM_trait_plot(subdat23,subdat23$SRL,"Specific root length") +
  geom_point(aes(y=quantile(SRL,.67),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = srl23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=1.5, vjust=-1.75, col="black") #add tukey labels

summary(rmf23.mod <- aov(RMF~trt, subdat23))
rmf23.tuk <- generateTukeyLabel(rmf23.mod, subdat23$RMF)
rmf.ca.23 <- CWM_trait_plot(subdat23,subdat23$RMF,"Root mass fraction") +
  geom_point(aes(y=quantile(RMF,.33),x="ir"), data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = rmf23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(lma23.mod <- aov(LMA~trt, subdat23))
lma23.tuk <- generateTukeyLabel(lma23.mod, subdat23$LMA)
lma.ca.23 <- CWM_trait_plot(subdat23,subdat23$LMA,"Leaf mass per area") +
  geom_point(aes(y=quantile(LMA,.67),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = lma23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") #add tukey labels

summary(seedmass23.mod <- aov(seed.mass~trt, subdat23))
seedmass23.tuk <- generateTukeyLabel(seedmass23.mod, subdat23$seed.mass)
seedmass.ca.23 <- CWM_trait_plot(subdat23,subdat23$seed.mass,"Seed mass") +
  #ylim(0,max(subdat23$seed.mass))+
  geom_point(aes(y=quantile(seed.mass,.67),x="dt"),data=traits.ca, col= "red", shape=18, size=3.5) +
  geom_text(data = seedmass23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=-1.75, col="black") 

# rootdiam.wy.23 <- CWM_trait_plot(cwm23.wy,cwm23.wy$rootdiam,"Root diameter")
#plot veg as Raoq instead (normalize FD/RoaQ)
subdatFD23$rootdiam <- normalize(subdatFD23$rootdiam) #(normalize FD/RoaQ)
summary(rd23.mod <- aov(rootdiam~trt, subdatFD23)) #summary model
rd23.tuk <- generateTukeyLabel(rd23.mod, subdatFD23$rootdiam) #tukey and labels
rootdiam.ca.23 <- CWM_trait_FD_plot(subdatFD23,subdatFD23$rootdiam,"Root diameter") +
  geom_point(aes(y=max(rootdiam),x="dt"),data=subdatFD23, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = rd23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=2.5, vjust=2, col="black") #add tukey labels

subdatFD23$full <- normalize(subdatFD23$full) #(normalize FD/RoaQ)
summary(full23.mod <- aov(full~trt, subdatFD23)) #summary model
full23.tuk <- generateTukeyLabel(full23.mod, subdatFD23$full) #tukey and labels
FD.ca.23 <- subdatFD23 %>% 
  ggplot(aes(trt,full)) +
  geom_violin(aes(fill=trt), width=1, trim=F) +
  geom_boxplot(width=.25) +
  #geom_jitter(width=.12, height = 0, size= 1, alpha=.3) +
  scale_fill_viridis_d(option="D", begin = .1, end = 1, alpha=.7) +
  theme_classic() + ylab("FD (RaoQ) of all traits") + xlab("") + 
  ylim(0,1)+
  #ylim(0,max(cwm_roaq23.ca$full+3)) +
  theme(legend.position  = "none") +
  geom_point(aes(y=max(full),x="fd"),data=subdatFD23, col= "red", shape=18, size=3.5) + #get target FD from annual highest possible (?)
  geom_text(data = full23.tuk, aes(x = trt, y = y_var, label = Letters), hjust=3, vjust=1, col="black") #add tukey labels


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

## figure for short report
#21
tiff("figures/cwm CA/alltargets_2021.tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.ca.21 + srl.ca.21 + rmf.ca.21) / (lma.ca.21 + seedmass.ca.21 + rootdiam.ca.21)) | (FD.ca.21)) +
  plot_layout(widths = c(2,1)) + 
  plot_annotation(title = 'CA 2021')
dev.off()
#22
tiff("figures/cwm CA/alltargets_2022.tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.ca.22 + srl.ca.22 + rmf.ca.22) / (lma.ca.22 + seedmass.ca.22 + rootdiam.ca.22)) | (FD.ca.22)) +
  plot_layout(widths = c(2,1)) + 
  plot_annotation(title = 'CA 2022')
dev.off()
#23
tiff("figures/cwm CA/alltargets_2023.tiff", res=400, height = 6,width =12, "in",compression = "lzw")
(((leafn.ca.23 + srl.ca.23 + rmf.ca.23) / (lma.ca.23 + seedmass.ca.23 + rootdiam.ca.23)) | (FD.ca.23)) +
  plot_layout(widths = c(2,1)) + 
  plot_annotation(title = 'CA 2023')
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