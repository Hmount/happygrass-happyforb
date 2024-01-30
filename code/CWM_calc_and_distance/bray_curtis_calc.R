#### Bray-curtis dissimilarity of community composition annually at both sites.

library(tidyverse)
library(vegan)


### CA
#get comp data
comp.ca <- read.csv("data/comp_ca.csv") #read in California comp data
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$Year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
comp.ca$trt <- tolower(comp.ca$trt) #make these lower to match cwm dataframe
comp.ca <- comp.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

# # run on all communities comparing to target community(same code)
# comms.ca <- comp.ca %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
# comms.ca <- comms.ca %>% column_to_rownames("trt.b.y") #  comp.ca[,c(18:56)] 
# comms.ca <- comms.ca[,c(16:54)]
# 
# # make comms.ca the long to add 4-letter codes
# codes <- data.frame(sixletter = colnames(comms.ca))
# codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
# colnames(comms.ca) <- codes$sppcodes
# colnames(comms.ca)[colnames(comms.ca) == "CACI"] <- "CAME" #change name of CACI to CAME

## pretreatment/ seeding probability 
# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
preds.ca <- read.csv("data/calgrass.allplot.assemblages.csv") #data
preds.ca$trt.b <- paste(preds.ca$trt, preds.ca$block)
preds.ca <- preds.ca %>% arrange(trt.b,block)
comms_p.ca <- labdsv::matrify(data.frame(preds.ca$trt.b,preds.ca$species,preds.ca$prob))
comms_p.ca <- comms_p.ca[,order(colnames(comms_p.ca))]

# #combine
# allcomms <- bind_rows(comms.ca,comms_p.ca)
# allcomms <- allcomms[,c(1:32)] #only OG 25 right now
# 
# #dissimilarity
# bcdist.ca <- vegdist.ca(as.matrix(allcomms), method="bray")
# bcdist.ca <-as.matrix(bcdist.ca)

##2021
# within each year subset the data
dist.ca21 <- comp.ca %>% filter(year=="2021")
dist.ca21 <- dist.ca21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
dist.ca21 <- dist.ca21 %>% column_to_rownames("trt.b.y")
dist.ca21 <- dist.ca21[,c(16:54)]
# make comms.ca the long to add 4-letter codes
codes <- data.frame(sixletter = colnames(dist.ca21))
codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
colnames(dist.ca21) <- codes$sppcodes
colnames(dist.ca21)[colnames(dist.ca21) == "CACI"] <- "CAME" #change name of CACI to CAME
dist.ca21 <- dist.ca21[ order(row.names(dist.ca21)), ]
# attach row of targets
comms_p.ca <- comms_p.ca[-c(99,205),] #remove "fd 50", not in comps data
all21 <- bind_rows(comms_p.ca, dist.ca21)
all21 <- all21[,c(1:32)] #only OG 25 right now
#run dist.ca or vegdist.ca
bcdist.camat <- vegdist(as.matrix(all21),method = "bray")
bcdist.camat <-as.matrix(bcdist.camat)[c(1:210),]
bcdist.camat <-as.matrix(bcdist.camat)[,-c(1:210)]
# save only pairwise between target 
bcdist.ca21<- data.frame(
  dist.ca=diag(as.matrix(bcdist.camat)),
  id=colnames(bcdist.camat))

##2022
# within each year subset the data
dist.ca22 <- comp.ca %>% filter(year=="2022")
dist.ca22 <- dist.ca22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
dist.ca22 <- dist.ca22 %>% column_to_rownames("trt.b.y")
dist.ca22 <- dist.ca22[,c(16:54)]
# make comms.ca the long to add 4-letter codes
#codes <- data.frame(sixletter = colnames(dist.ca22))
#codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
colnames(dist.ca22) <- codes$sppcodes
colnames(dist.ca22)[colnames(dist.ca22) == "CACI"] <- "CAME" #change name of CACI to CAME
dist.ca22 <- dist.ca22[ order(row.names(dist.ca22)), ]
# attach row of targets
#comms_p.ca <- comms_p.ca[-c(99,205),] #remove "ir 50" and "fd 50", not in comps data
all22 <- bind_rows(comms_p.ca, dist.ca22)
all22 <- all22[,c(1:32)] #only OG 25 right now
#run dist.ca or vegdist.ca
bcdist.camat <- vegdist(as.matrix(all22),method = "bray")
bcdist.camat <-as.matrix(bcdist.camat)[c(1:210),]
bcdist.camat <-as.matrix(bcdist.camat)[,-c(1:210)]
# save only pairwise between target 
bcdist.ca22<- data.frame(
  dist.ca=diag(as.matrix(bcdist.camat)),
  id=colnames(bcdist.camat))

##2023
# within each year subset the data
dist.ca23 <- comp.ca %>% filter(year=="2023")
dist.ca23 <- dist.ca23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
dist.ca23 <- dist.ca23 %>% column_to_rownames("trt.b.y")
dist.ca23 <- dist.ca23[,c(16:54)]
# make comms.ca the long to add 4-letter codes
#codes <- data.frame(sixletter = colnames(dist.ca23))
#codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
colnames(dist.ca23) <- codes$sppcodes
colnames(dist.ca23)[colnames(dist.ca23) == "CACI"] <- "CAME" #change name of CACI to CAME
dist.ca23 <- dist.ca23[ order(row.names(dist.ca23)), ]
# attach row of targets
#comms_p.ca <- comms_p.ca[-c(99,205),] #remove "ir 50" and "fd 50", not in comps data
all23 <- bind_rows(comms_p.ca, dist.ca23)
all23 <- all23[,c(1:32)] #only OG 25 right now
#run dist.ca or vegdist.ca
bcdist.camat <- vegdist(as.matrix(all23),method = "bray")
bcdist.camat <-as.matrix(bcdist.camat)[c(1:210),]
bcdist.camat <-as.matrix(bcdist.camat)[,-c(1:210)]
# save only pairwise between target 
bcdist.ca23<- data.frame(
  dist.ca=diag(as.matrix(bcdist.camat)),
  id=colnames(bcdist.camat))

## combine dissimilarity dataframe
bcdist.ca <- bind_rows(bcdist.ca21,bcdist.ca22,bcdist.ca23)

#save
write.csv(bcdist.ca, "data/bc_dissimilarity_ca.csv" )


### WY
## load in composition data, clean and modify columns as usual
comp.wy <- read.csv("data/raw_cover/hpg_total.csv") #Wyoming species comp data 
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
######
##FIX PROPORTIONS AND PROPGATE THROUGHOUT
#### USE PLOT LEVEL:
comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data


## pretreatment/ seeding probability 
# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
preds.wy <- read.csv("data/allplot.assemblages.csv") #data
preds.wy <- preds.wy %>% arrange(trt,block)
comms_p.wy <- labdsv::matrify(data.frame(preds.wy$trt,preds.wy$species,preds.wy$prob))
comms_p.wy <- comms_p.wy[,order(colnames(comms_p.wy))]

# #combine
# allcomms <- bind_rows(comms.ca,comms_p.ca)
# allcomms <- allcomms[,c(1:32)] #only OG 25 right now
# 
# #dissimilarity
# bcdist.ca <- vegdist.ca(as.matrix(allcomms), method="bray")
# bcdist.ca <-as.matrix(bcdist.ca)

##2021
# within each year subset the data
dist.wy21 <- comp.wy.wide %>% filter(year=="2021")
dist.wy21 <- dist.wy21 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
dist.wy21 <- dist.wy21 %>% column_to_rownames("trt.b.sub.y")
dist.wy21 <- dist.wy21[,c(7:62)]
# make comms.wy the long to add 4-letter codes
dist.wy21 <- dist.wy21[ order(row.names(dist.wy21)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all21 <- bind_rows(comms_p.wy, dist.wy21)
all21 <- all21[,c(1:25)] #only OG 25 right now
all21[is.na(all21)] <- 0
#run dist.wy or vegdist.wy
bcdist.wymat <- vegdist(as.matrix(all21),method = "bray")
bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
#north first
filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
  select(matches("*\\.n\\.*"))
# save only pairwise between target 
bcdist.wy21.n<- data.frame(
  dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#south next
filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
  select(matches("*\\.s\\.*"))
# save only pairwise between target 
bcdist.wy21.s<- data.frame(
  dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#combine
bcdist.wy21 <- bind_rows(bcdist.wy21.n,bcdist.wy21.s)

##2022
# within each year subset the data
dist.wy22 <- comp.wy.wide %>% filter(year=="2022")
dist.wy22 <- dist.wy22 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
dist.wy22 <- dist.wy22 %>% column_to_rownames("trt.b.sub.y")
dist.wy22 <- dist.wy22[,c(7:62)]
# make comms.wy the long to add 4-letter codes
dist.wy22 <- dist.wy22[ order(row.names(dist.wy22)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all22 <- bind_rows(comms_p.wy, dist.wy22)
all22 <- all22[,c(1:25)] #only OG 25 right now
all22[is.na(all22)] <- 0
#run dist.wy or vegdist.wy
bcdist.wymat <- vegdist(as.matrix(all22),method = "bray")
bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
#north first
filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
  select(matches("*\\.n\\.*"))
# save only pairwise between target 
bcdist.wy22.n<- data.frame(
  dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#south next
filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
  select(matches("*\\.s\\.*"))
# save only pairwise between target 
bcdist.wy22.s<- data.frame(
  dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#combine
bcdist.wy22 <- bind_rows(bcdist.wy22.n,bcdist.wy22.s)

##2023
# within each year subset the data
dist.wy23 <- comp.wy.wide %>% filter(year=="2023")
dist.wy23 <- dist.wy23 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
dist.wy23 <- dist.wy23 %>% column_to_rownames("trt.b.sub.y")
dist.wy23 <- dist.wy23[,c(7:62)]
# make comms.wy the long to add 4-letter codes
dist.wy23 <- dist.wy23[ order(row.names(dist.wy23)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all23 <- bind_rows(comms_p.wy, dist.wy23)
all23 <- all23[,c(1:25)] #only OG 25 right now
all23[is.na(all23)] <- 0
#run dist.wy or vegdist.wy
bcdist.wymat <- vegdist(as.matrix(all23),method = "bray")
bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
#north first
filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
  select(matches("*\\.n\\.*"))
# save only pairwise between target 
bcdist.wy23.n<- data.frame(
  dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#south next
filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
  select(matches("*\\.s\\.*"))
# save only pairwise between target 
bcdist.wy23.s<- data.frame(
  dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#combine
bcdist.wy23 <- bind_rows(bcdist.wy23.n,bcdist.wy23.s)

## combine dissimilarity dataframe
bcdist.wy <- bind_rows(bcdist.wy21,bcdist.wy22,bcdist.wy23)

#save
write.csv(bcdist.wy, "data/bc_dissimilarity_wy.csv" )

#### figures 
#CA
x <- read.csv("data/comp_ca.csv")
x$trt <- tolower(x$trt)
x <- x %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
droughtca <- x %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=F)
droughtca <- droughtca %>% select(trt.b.y,water)
bcdist.ca <- merge(bcdist.ca,droughtca, by.y="trt.b.y",by.x="id", all.x=T)
bcdist.ca <- separate(bcdist.ca, id, into = c("trt", "block", "year"), sep = "\\.")
bcdist.ca$year <- as.factor(bcdist.ca$year)
bcdist.ca$block <- as.factor(bcdist.ca$block)
bcdist.ca$trt <- as.factor(bcdist.ca$trt)
bcdist.ca$water <- as.factor(bcdist.ca$water)
bcdist.ca <- bcdist.ca %>% mutate(drought=as.factor(ifelse(water=="0.5","drt","cntl")))
droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

summary(aov(dist.ca~trt*drought*year, bcdist.ca))
ggplot(bcdist.ca, aes(y=dist.ca, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values=droughtcolsca)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()

#WY
x <- read.csv("data/comp_wy.csv")
x <- x %>% filter(year!="2020")
droughtwy <- x %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=F)
droughtwy <- droughtwy %>% select(trt.b.sub.y,drought)
bcdist.wy <- merge(bcdist.wy,droughtwy, by.y="trt.b.sub.y",by.x="id", all.x=T)
bcdist.wy <- separate(bcdist.wy, id, into = c("trt", "block", "subplot","year"), sep = "\\.")
bcdist.wy$year <- as.factor(bcdist.wy$year)
bcdist.wy$block <- as.factor(bcdist.wy$block)
bcdist.wy$trt <- as.factor(bcdist.wy$trt)
bcdist.wy <- bcdist.wy %>% mutate(drought=as.factor(ifelse(drought=="0","cntl","drt")))
#bcdist.wy$drought <- as.factor(bcdist.wy$drought)
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

summary(aov(dist.wy~trt*drought*year, bcdist.wy))
ggplot(bcdist.wy, aes(y=dist.wy, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values=droughtcolswy)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()
