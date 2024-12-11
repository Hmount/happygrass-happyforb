#### Relative cover of all species in each community by functional group

library(tidyverse)

# load community comp data
compwy <- read.csv("data/comp_wy_plot.csv")
compca <- read.csv("data/comp_ca.csv")

### summarize the total per functional group (WY) 
### no breaking it down this way currently
# compwy2 <- unite(compwy, trt.b.y, c("trt", "block", "year"), remove = F) #make plot ID column
# subdat <- compwy2 %>%
#   group_by(trt.b.y, native, graminoid) %>%  # Group by native status and graminoid status
#   mutate(groupcov = sum(cover, na.rm = TRUE)) # Sum cover within each group
# subdat2 <- subdat %>%
#   filter(year!="2020") %>%
#   group_by(year, native, graminoid) %>%  # Group by native status and graminoid status
#   summarize(groupmeancov = mean(groupcov, na.rm = TRUE))
# subdat2$groupcol <- c("I0","I1","N0","N1","na","I0","I1","N0","N1","na","I0","I1","N0","N1","na")
# subdat3 <- subdat2 %>% filter(groupcol == "I1" | groupcol == "N1"| groupcol == "N0")
# ggplot(subdat3, aes(y=groupmeancov, x=year, fill=groupcol))+
#   geom_col()

# group within year and precip treatment, find plot-level averages, prep data 
tttw <- compwy %>% group_by(year,drought) %>%
  summarise(nativemean = mean(nativecov.plot),
            litmean = mean(lit.plotmean, na.rm=T),
            bgmean = mean(bg.plotmean, na.rm=T),
            brtemean = mean(BRTE))
tttw <- tttw %>% filter(year!="2020") #remove pre-treatment (2020)
#pivot for plotting
tttwwide <- tttw %>% pivot_longer(cols= c(nativemean,litmean,bgmean,brtemean), names_to = "covtype", values_to = "coverprop")
tttwwide$drought <- as.factor(tttwwide$drought) #make drought a factor
tttwwide <- tttwwide %>% mutate(coverper=coverprop*100)
tttwwide$drought2 <- factor(tttwwide$drought, levels = c("1", "0"))
#make a nice plot
wycovplot <- ggplot(tttwwide, aes(x=drought2,y=coverper, fill=covtype))+
  geom_col()+
  scale_fill_manual(values=c("grey60","salmon", "gold","darkgreen"),
                    labels=c("Bare ground","Invasive grass","Litter","Native cover"))+
  facet_wrap(~year)+
  scale_x_discrete(labels = c("Reduction", "Ambient"))+
  labs(y="Absolute cover (%)", 
       x="Precipitation treatment", 
       fill="Cover type",
       title = "Wyoming")+
  theme_classic2()+
  theme(axis.text.x= element_text(angle=15, hjust=.5))


### can we do CA?
comp.ca <- read.csv("data/raw_cover/Species_Composition_allyears.csv") #read in California comp data
comp.ca$bgunknown <- 1-comp.ca$native.cover-comp.ca$inv.grass.cov
comp.ca$bgunknown[comp.ca$bgunknown <= 0] <- 0 #make negatives into 0 
# group within year and precip treatment, find plot-level averages, prep data 
tttc <- comp.ca %>% group_by(Year,water) %>%
  summarise(nativemean = mean(native.cover),
            PHCImean = mean(PHACIC, na.rm=T),
            AMMEmean = mean(AMSMEN, na.rm=T),
            a.bgmean = mean(bgunknown, na.rm=T),
            b.invmean = mean(inv.grass.cov))
tttc <- tttc %>% mutate(c.weedyforbs = PHCImean+AMMEmean)
tttc <- tttc %>% mutate(d.nativeadjusted = nativemean-c.weedyforbs)
#pivot for plotting
tttcwide <- tttc %>% pivot_longer(cols= c(d.nativeadjusted,c.weedyforbs,a.bgmean,b.invmean), names_to = "covtype", values_to = "coverprop")
tttcwide$water <- as.factor(tttcwide$water) #make water a factor
tttcwide <- tttcwide %>% mutate(coverper=coverprop*100)
#make a nice plot
cacovplot <-ggplot(tttcwide, aes(x=water,y=coverper, fill=covtype))+
  geom_col()+
  scale_fill_manual(values=c("grey60","salmon","gold","darkgreen"),
                    labels=c("Bare ground","Invasive grass","Weedy forbs","Native cover"))+
  facet_wrap(~Year)+
  scale_x_discrete(labels = c("Reduction", "Addition"))+
  labs(y="Absolute cover (%)", 
       x="Precipitation treatment", 
       fill="Cover type",
       title = "California")+
  theme_classic2()+
  theme(axis.text.x= element_text(angle=15, hjust=.5))


library(ggpubr)
#export for short report
tiff("figures/covplots.tiff", res=400, height = 6,width =6, "in",compression = "lzw")
ggarrange(wycovplot,cacovplot, nrow=2)
dev.off()