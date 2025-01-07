#### Bray-curtis dissimilarity of community composition annually at both sites.

#### Packages
library(tidyverse)
library(vegan)

#### Calculate B-C dissimilarity for every community every year
####
#### CA
####
## load in composition data, clean and modify columns as usual
comp.ca <- read.csv("data/comp_ca.csv") #read in California comp data
comp.ca$fesper.seeded <- as.factor(comp.ca$fesper.seeded)
comp.ca$fesper.present <- as.factor(comp.ca$fesper.present)
comp.ca$water <- as.factor(comp.ca$water)
comp.ca$year<-as.factor(comp.ca$year)
comp.ca$trt <- as.factor(comp.ca$trt)
comp.ca$block <- as.factor(comp.ca$block)
#comp.ca$trt <- tolower(comp.ca$trt) #make these lower to match cwm dataframe
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
dist.ca21 <- dist.ca21[,c(15:53)]
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
## Should remove rows with no species since we can't compare a community to no community
## but these comms have spp (just not the native comm, so 0 comparison is valid enough
## for my purposes at the moment (2/2/24)
#test <- all21 %>% filter(rowSums(all21)=="0") 
 # if do this, also remove the matching target rows
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
dist.ca22 <- dist.ca22[,c(15:53)]
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
all22 <- all22 %>% filter(rowSums(all22) != "0")
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
dist.ca23 <- dist.ca23[,c(15:53)]
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
all23 <- all23 %>% filter(rowSums(all23) != "0")
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
write.csv(bcdist.ca, "data/bc_dissimilarity_ca.csv", row.names = F)

####
#### WY
####
## load in composition data, clean and modify columns as usual
comp.wy.plot <- read.csv("data/comp_wy_plot.csv")
comp.wy.plot <- comp.wy.plot %>% filter(year != "2020") #remove 2020
comp.wy.plot$drought <- as.factor(comp.wy.plot$drought)
comp.wy.plot$year<-as.factor(comp.wy.plot$year)
comp.wy.plot$trt <- as.factor(comp.wy.plot$trt)
comp.wy.plot$block <- as.factor(comp.wy.plot$block)
#comp.wy.plot <- comp.wy.plot %>% unite(plot, c(block, trt), sep = ".", remove=F) # make unique plot variable

## pre-treatment/ seeding probability
# Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
preds.wy <- read.csv("data/allplot.assemblages.csv") #data
preds.wy <- preds.wy %>% arrange(trt,block)
comms_p.wy <- labdsv::matrify(data.frame(preds.wy$trt,preds.wy$species,preds.wy$prob))
comms_p.wy <- comms_p.wy[,order(colnames(comms_p.wy))]

##2021
# within each year subset the data
dist.wy21 <- comp.wy.plot %>% filter(year=="2021")
dist.wy21 <- dist.wy21 %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T)
dist.wy21 <- dist.wy21 %>% column_to_rownames("trt.b.y")
dist.wy21 <- dist.wy21[,c(3:58)]
# order comms
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
# save only pairwise between target
bcdist.wy21<- data.frame(
  dist.wy=diag(as.matrix(bcdist.wymat)),
  id=colnames(bcdist.wymat))

##2022
# within each year subset the data
dist.wy22 <- comp.wy.plot %>% filter(year=="2022")
dist.wy22 <- dist.wy22 %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T)
dist.wy22 <- dist.wy22 %>% column_to_rownames("trt.b.y")
dist.wy22 <- dist.wy22[,c(3:58)]
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
# save only pairwise between target
bcdist.wy22<- data.frame(
  dist.wy=diag(as.matrix(bcdist.wymat)),
  id=colnames(bcdist.wymat))

##2023
# within each year subset the data
dist.wy23 <- comp.wy.plot %>% filter(year=="2023")
dist.wy23 <- dist.wy23 %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T)
dist.wy23 <- dist.wy23 %>% column_to_rownames("trt.b.y")
dist.wy23 <- dist.wy23[,c(3:58)]
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
# save only pairwise between target
bcdist.wy23<- data.frame(
  dist.wy=diag(as.matrix(bcdist.wymat)),
  id=colnames(bcdist.wymat))

## combine dissimilarity dataframe
bcdist.wy <- bind_rows(bcdist.wy21,bcdist.wy22,bcdist.wy23)

#save
write.csv(bcdist.wy, "data/bc_dissimilarity_wy.csv", row.names = F)

##2021
# within each year subset the data
test <- comp.wy.plot %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=T)
test <- test %>% column_to_rownames("trt.b.y")
test <- test[,c(3:58)]
# order comms
test <- test[ order(row.names(dist.wy21)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all21 <- bind_rows(comms_p.wy, dist.wy21)
all21 <- all21[,c(1:25)] #only OG 25 right now
test[is.na(test)] <- 0
#run dist.wy or vegdist.wy
testmat <- vegdist(as.matrix(test),method = "bray")
test2<-data.frame(as.matrix(testmat))
bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
# save only pairwise between target
bcdist.wy21<- data.frame(
  dist.wy=diag(as.matrix(bcdist.wymat)),
  id=colnames(bcdist.wymat))


# ### WY (PREVIOUS WITH INCORRECT % DATA, 
######### could be used again to get subplot dissimilarity if we end up wanting it)
# ## load in composition data, clean and modify columns as usual
# comp.wy <- read.csv("data/raw_cover/hpg_total.csv") #Wyoming species comp data 
# comp.wy <- comp.wy %>% filter(year != "2020") #keep only 2023 data 
# comp.wy$subplot <- as.factor(comp.wy$subplot)
# comp.wy$drought <- as.factor(comp.wy$drought)
# comp.wy$year<-as.factor(comp.wy$year)
# comp.wy$trt <- as.factor(comp.wy$trt)
# comp.wy$block <- as.factor(comp.wy$block)
# comp.wy <- comp.wy %>% unite(plot, c(block, trt, subplot), sep = ".", remove=F) # make unique plot variable
# fornativecover <- comp.wy %>% filter(species!="BG"&
#                                        species!="Litter"&
#                                        native == "N") %>% #only native live veg
#   group_by(year,block,trt,subplot) %>% 
#   summarize(nativecov = sum(cover, na.rm=T)) #summarize total live veg per subplot
# comp.wy <- merge(comp.wy,fornativecover, all.x = T)
# comp.wy.wide <- comp.wy %>% select(-c("prob","native","graminoid")) #columns to drop 
# comp.wy.wide <- comp.wy.wide %>% pivot_wider(id_cols = c("year","block","trt","subplot","drought",
#                                                          "nativecov","BG", "Litter","plot","sub.tveg"), 
#                                              names_from = "species", 
#                                              values_from = "cover")
# 
# comp.wy.wide$nativecov <- comp.wy.wide$nativecov/100  # make native live veg % a proportion to match CA data
# 
# 
# ## pretreatment/ seeding probability 
# # Calculating community weighted means for the seeded communities. Needed for determine proximity to our objective.
# preds.wy <- read.csv("data/allplot.assemblages.csv") #data
# preds.wy <- preds.wy %>% arrange(trt,block)
# comms_p.wy <- labdsv::matrify(data.frame(preds.wy$trt,preds.wy$species,preds.wy$prob))
# comms_p.wy <- comms_p.wy[,order(colnames(comms_p.wy))]
# 
# # #combine
# # allcomms <- bind_rows(comms.ca,comms_p.ca)
# # allcomms <- allcomms[,c(1:32)] #only OG 25 right now
# # 
# # #dissimilarity
# # bcdist.ca <- vegdist.ca(as.matrix(allcomms), method="bray")
# # bcdist.ca <-as.matrix(bcdist.ca)
# 
# ##2021
# # within each year subset the data
# dist.wy21 <- comp.wy.wide %>% filter(year=="2021")
# dist.wy21 <- dist.wy21 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
# dist.wy21 <- dist.wy21 %>% column_to_rownames("trt.b.sub.y")
# dist.wy21 <- dist.wy21[,c(7:62)]
# # make comms.wy the long to add 4-letter codes
# dist.wy21 <- dist.wy21[ order(row.names(dist.wy21)), ]
# # attach row of targets
# #comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
# all21 <- bind_rows(comms_p.wy, dist.wy21)
# all21 <- all21[,c(1:25)] #only OG 25 right now
# all21[is.na(all21)] <- 0
# #run dist.wy or vegdist.wy
# bcdist.wymat <- vegdist(as.matrix(all21),method = "bray")
# bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
# bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
# #north first
# filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
#   select(matches("*\\.n\\.*"))
# # save only pairwise between target 
# bcdist.wy21.n<- data.frame(
#   dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
#   id=colnames(filtered_dissimilarity_df))
# #south next
# filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
#   select(matches("*\\.s\\.*"))
# # save only pairwise between target 
# bcdist.wy21.s<- data.frame(
#   dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
#   id=colnames(filtered_dissimilarity_df))
# #combine
# bcdist.wy21 <- bind_rows(bcdist.wy21.n,bcdist.wy21.s)
# 
# ##2022
# # within each year subset the data
# dist.wy22 <- comp.wy.wide %>% filter(year=="2022")
# dist.wy22 <- dist.wy22 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
# dist.wy22 <- dist.wy22 %>% column_to_rownames("trt.b.sub.y")
# dist.wy22 <- dist.wy22[,c(7:62)]
# # make comms.wy the long to add 4-letter codes
# dist.wy22 <- dist.wy22[ order(row.names(dist.wy22)), ]
# # attach row of targets
# #comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
# all22 <- bind_rows(comms_p.wy, dist.wy22)
# all22 <- all22[,c(1:25)] #only OG 25 right now
# all22[is.na(all22)] <- 0
# #run dist.wy or vegdist.wy
# bcdist.wymat <- vegdist(as.matrix(all22),method = "bray")
# bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
# bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
# #north first
# filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
#   select(matches("*\\.n\\.*"))
# # save only pairwise between target 
# bcdist.wy22.n<- data.frame(
#   dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
#   id=colnames(filtered_dissimilarity_df))
# #south next
# filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
#   select(matches("*\\.s\\.*"))
# # save only pairwise between target 
# bcdist.wy22.s<- data.frame(
#   dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
#   id=colnames(filtered_dissimilarity_df))
# #combine
# bcdist.wy22 <- bind_rows(bcdist.wy22.n,bcdist.wy22.s)
# 
# ##2023
# # within each year subset the data
# dist.wy23 <- comp.wy.wide %>% filter(year=="2023")
# dist.wy23 <- dist.wy23 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
# dist.wy23 <- dist.wy23 %>% column_to_rownames("trt.b.sub.y")
# dist.wy23 <- dist.wy23[,c(7:62)]
# # make comms.wy the long to add 4-letter codes
# dist.wy23 <- dist.wy23[ order(row.names(dist.wy23)), ]
# # attach row of targets
# #comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
# all23 <- bind_rows(comms_p.wy, dist.wy23)
# all23 <- all23[,c(1:25)] #only OG 25 right now
# all23[is.na(all23)] <- 0
# #run dist.wy or vegdist.wy
# bcdist.wymat <- vegdist(as.matrix(all23),method = "bray")
# bcdist.wymat <-as.matrix(bcdist.wymat)[c(1:256),]
# bcdist.wymat <-as.matrix(bcdist.wymat)[,-c(1:256)]
# #north first
# filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
#   select(matches("*\\.n\\.*"))
# # save only pairwise between target 
# bcdist.wy23.n<- data.frame(
#   dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
#   id=colnames(filtered_dissimilarity_df))
# #south next
# filtered_dissimilarity_df <- data.frame(bcdist.wymat) %>%
#   select(matches("*\\.s\\.*"))
# # save only pairwise between target 
# bcdist.wy23.s<- data.frame(
#   dist.wy=diag(as.matrix(filtered_dissimilarity_df)),
#   id=colnames(filtered_dissimilarity_df))
# #combine
# bcdist.wy23 <- bind_rows(bcdist.wy23.n,bcdist.wy23.s)
# 
# ## combine dissimilarity dataframe
# bcdist.wy <- bind_rows(bcdist.wy21,bcdist.wy22,bcdist.wy23)
# 
# #save
# write.csv(bcdist.wy, "data/bc_dissimilarity_wy.csv" )


#### models and figures ####
####
#### CA
####
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

summary(bcca<-aov(dist.ca~trt*drought*year, bcdist.ca))
tukbcca <- TukeyHSD(bcca)

letterstest <- data.frame(multcompView::multcompLetters4(bcca,tukbcca)$'trt:drought:year'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\2", rownames(letterstest)))
letterstest$year <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\3", rownames(letterstest)))
test <- bcdist.ca %>% group_by(year, drought, trt) %>% summarise(yposition = quantile(dist.ca,.8))
test <- merge(letterstest,test, by = c("year", "drought", "trt"))
testbc <- merge(test,bcdist.ca, by = c("year", "drought", "trt"), all=T)

ggplot(testbc, aes(y=dist.ca, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  scale_fill_manual(values=droughtcolsca)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()

#export for short report
tiff("figures/cwm ca/dissim_ca.tiff", res=400, height = 4,width =8.5, "in",compression = "lzw")
ggplot(testbc, aes(y=dist.ca, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  scale_fill_manual(values=droughtcolsca)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()
dev.off()



#WY
x <- read.csv("data/comp_wy_plot.csv")
x <- x %>% filter(year!="2020")
droughtwy <- x %>% unite(trt.b.y, c(trt, block,year), sep = ".", remove=F)
droughtwy <- droughtwy %>% select(trt.b.y,drought)
bcdist.wy <- merge(bcdist.wy,droughtwy, by.y="trt.b.y",by.x="id", all.x=T)
bcdist.wy <- separate(bcdist.wy, id, into = c("trt", "block","year"), sep = "\\.")
bcdist.wy$year <- as.factor(bcdist.wy$year)
bcdist.wy$block <- as.factor(bcdist.wy$block)
bcdist.wy$trt <- as.factor(bcdist.wy$trt)
bcdist.wy <- bcdist.wy %>% mutate(drought=as.factor(ifelse(drought=="0","cntl","drt")))
#bcdist.wy$drought <- as.factor(bcdist.wy$drought)
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

summary(bcwy<-aov(dist.wy~trt*drought*year, bcdist.wy))
tukbcwy <- TukeyHSD(bcwy)

letterstest <- data.frame(multcompView::multcompLetters4(bcwy,tukbcwy)$'trt:drought:year'['Letters'])
letterstest$trt <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\1", rownames(letterstest)))
letterstest$drought <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\2", rownames(letterstest)))
letterstest$year <- as.factor(sub("([a-z]+):([a-z]+):(\\d+)$", "\\3", rownames(letterstest)))
test <- bcdist.wy %>% group_by(year, drought, trt) %>% summarise(yposition = quantile(dist.wy,.8))
test <- merge(letterstest,test, by = c("year", "drought", "trt"))
testbc <- merge(test,bcdist.wy, by = c("year", "drought", "trt"), all=T)

ggplot(testbc, aes(y=dist.wy, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  scale_fill_manual(values=droughtcolswy)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()

#export for short report
tiff("figures/cwm wy/dissim_wy.tiff", res=400, height = 4,width =8.5, "in",compression = "lzw")
ggplot(testbc, aes(y=dist.wy, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  geom_text(aes(y=yposition,label = Letters), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            #angle = 15,
            size=3) +
  scale_fill_manual(values=droughtcolswy)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()
dev.off()



###### TESTING/CHECKING B-C and FD
FDdat <- read.csv("data/cwm_raoq_wy(plot).csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)

bcdist.wy <- read.csv("data/bc_dissimilarity_wy.csv") #read in comp data
bcdist.wy <- separate(bcdist.wy, id, into = c("trt", "block", "year"), sep = "\\.")
bcdist.wy <- bcdist.wy %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

test.wy <- merge(bcdist.wy,FDdat, by= c("trt", "block", "year"))
subvalid <- comp.wy %>% group_by(block,trt,year) %>%
  mutate(propnative = nativecov.plot/totcov.plot*100) %>% 
  filter(propnative < 80)
subvalid <- subvalid %>%  unite(trt.b.y, c(trt, block, year), sep = ".") 
test.wy <- test.wy %>%  unite(trt.b.y, c(trt, block, year), sep = ".")
testvalid <- test.wy %>% filter(!(trt.b.y %in% subvalid$trt.b.y))
testvalid <- separate(testvalid, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")

summary(lm(full~dist.wy*year*trt*drought, testvalid)) #signifigany and explains a lot
ggplot(testvalid, aes(x=dist.wy,y=full,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year, scales="free")
# more dissimilar in composition lead to lower functional diversities
# our targets were highly FD, and real/stable communities were less FD?

wydist2 <- read.csv("data/cwm_maxdistances_wy(plot).csv")

wydist2 <- separate(wydist2, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
test <- merge(bcdist.wy,wydist2, by= c("trt", "block", "year"), all.x=T)

subvalid <- comp.wy %>% group_by(block,trt,year) %>%
  mutate(propnative = nativecov.plot/totcov.plot*100) %>% 
  filter(propnative < 80)
subvalid <- subvalid %>%  unite(trt.b.y, c(trt, block, year), sep = ".") 
test <- test %>%  unite(trt.b.y, c(trt, block, year), sep = ".")
testvalid <- test %>% filter(!(trt.b.y %in% subvalid$trt.b.y))
testvalid <- separate(testvalid, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
subdrought <- FDdat %>% select(c(year,trt,block,drought))
testvalid <- merge(testvalid, subdrought, all.x=T)

#quick and dirty make this plot, make dist column that is dist to target
testvalid.wy <- testvalid %>% mutate(dist_target = ifelse(trt=="dt", distdt,
                                                          ifelse(trt=="ir", distir,
                                                                 ifelse(trt=="fd", distfd,distr))))

summary(lm(dist_target~dist.wy*trt*year*drought, testvalid)) #signifigany and explains a lot
ggplot(testvalid, aes(x=dist.wy,y=dist_target,col=trt))+
  geom_point()+
  geom_smooth(method="lm")#+
#facet_grid(year~drought, scales="free")#+
#facet_wrap(~year, scales = "free")
# mix of relationships, but positive slopes in FD and IR21+22 are expected relationships
# DT is oddly static and then negative. Random has no relationship (matches CA)

#export for short report
tiff("figures/cwm wy/dist_dis_wy.tiff", res=400, height = 4,width =6, "in",compression = "lzw")
ggplot(testvalid, aes(x=dist.wy,y=dist_target,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  #facet_grid(year~drought, scales="free")#+
  labs(x="Compositional dissimilarity (bray-curtis)", 
       y="Distance to functonal targets (Euclidean)",
       col="seeding trt")+
  theme_classic()
#facet_wrap(~year, scales = "free")
dev.off()


### CA
#what about with FD
FDdat <- read.csv("data/cwm_raoq_ca.csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)
FDdat$subplot <- as.factor(FDdat$water)

bcdist.ca <- read.csv("data/bc_dissimilarity_ca.csv") #read in comp data
bcdist.ca <- separate(bcdist.ca, id, into = c("trt", "block", "year"), sep = "\\.")
bcdist.ca <- bcdist.ca %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df

test.ca <- merge(bcdist.ca,FDdat, by= c("trt", "block", "year"))
subvalid <- comp.ca %>% group_by(block,trt,year) %>%
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>% 
  filter(propnative < 80)
subvalid <- subvalid %>%  unite(trt.b.y, c(trt, block, year), sep = ".") 
test.ca <- test.ca %>%  unite(trt.b.y, c(trt, block, year), sep = ".")
testvalid <- test.ca %>% filter(!(trt.b.y %in% subvalid$trt.b.y))
testvalid <- separate(testvalid, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")

summary(lm(full~dist.ca*year*trt*water, testvalid)) #signifigant but explains little
ggplot(testvalid, aes(x=dist.ca,y=full,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year)
# more dissimilar compositions lead to lower functional diversities 
# our targets were highly FD, and going away from that homogenized

## with distance and dissimilrity? 
cadist2 <- read.csv("data/cwm_maxdistances_ca.csv")

cadist2 <- separate(cadist2, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")
test <- merge(bcdist.ca,cadist2, by= c("trt", "block", "year"))
subvalid <- comp.ca %>% group_by(block,trt,year) %>%
  mutate(propnative = native.cover/(native.cover+inv.grass.cov)*100) %>% 
  filter(propnative < 80)
subvalid <- subvalid %>%  unite(trt.b.y, c(trt, block, year), sep = ".") 
test <- test %>%  unite(trt.b.y, c(trt, block, year), sep = ".")
testvalid <- test %>% filter(!(trt.b.y %in% subvalid$trt.b.y))
testvalid <- separate(testvalid, trt.b.y, into = c("trt", "block", "year"), sep = "\\.")

#quick and dirty make this plot, make dist column that is dist to target
testvalid.ca <- testvalid %>% mutate(dist_target = ifelse(trt=="dt", distdt,
                                                          ifelse(trt=="ir", distir,
                                                                 ifelse(trt=="fd", distfd,distr))))

summary(lm(dist_target~dist.ca*year*trt, testvalid.ca)) #significant but explains little
ggplot(testvalid.ca, aes(x=dist.ca,y=dist_target,col=trt))+
  geom_point()+
  geom_smooth(method="lm")#+
#facet_wrap(~year)
# more similarity to composition = more similarity to target (except random)
# this is the expected result

#export for short report
tiff("figures/cwm ca/dist_dis_ca.tiff", res=400, height = 4,width =8.5, "in",compression = "lzw")
ggplot(testvalid, aes(x=dist.ca,y=dist,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  #facet_grid(year~drought, scales="free")#+
  labs(x="Compositional dissimilarity (bray-curtis)", 
       y="Distance to functonal targets (Euclidean)",
       col="seeding trt")+
  facet_wrap(~year, scales = "free")
dev.off()

####together for 2/5/24
library(tidyverse)
library(patchwork)
testvalid.ca$trt <- factor(testvalid.ca$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
pdisdist.ca <- ggplot(testvalid.ca, aes(x=dist.ca,y=dist_target,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(x=" ", 
       y=" ",
       col="seeding trt")+
  scale_color_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7) +
  theme_classic()
testvalid.wy$trt <- factor(testvalid.wy$trt, levels = c('ir','dt','fd','rand'),ordered = TRUE) #order trt levels
pdisdist.wy <- ggplot(testvalid.wy, aes(x=dist.wy,y=dist_target,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  scale_color_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7) +
  #scale_color_viridis_d(option = "D", begin = 1, end = 0.1, alpha = 0.7) +
  labs(x=" ", 
       y=" ",
       col="seeding trt")+
  theme_classic()
library(ggpubr)
annotate_figure(ggarrange(pdisdist.ca,pdisdist.wy, common.legend = T),
                bottom ="Compositional dissimilarity (bray-curtis)",
                left = "Distance to functonal targets (Euclidean)")


# #export for short report
# tiff("figures/cwm wy/dissim_wy.tiff", res=400, height = 4,width =8.5, "in",compression = "lzw")
# ggplot(testbc, aes(y=dist., x=trt, fill=drought))+
#   geom_boxplot()+
#   #geom_smooth(method="lm")+
#   geom_text(aes(y=yposition,label = Letters), 
#             position = position_dodge(width = 0.9), 
#             vjust = -0.5,
#             #angle = 15,
#             size=3) +
#   scale_fill_manual(values=droughtcolswy)+
#   facet_wrap(~year, scales="fixed")+
#   labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
#   theme_classic()
# dev.off()

