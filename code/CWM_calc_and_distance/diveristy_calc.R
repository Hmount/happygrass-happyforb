#### Diveristy of community composition annually at both sites.

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
# preds.ca <- read.csv("data/calgrass.allplot.assemblages.csv") #data
# preds.ca$trt.b <- paste(preds.ca$trt, preds.ca$block)
# preds.ca <- preds.ca %>% arrange(trt.b,block)
# comms_p.ca <- labdsv::matrify(data.frame(preds.ca$trt.b,preds.ca$species,preds.ca$prob))
# comms_p.ca <- comms_p.ca[,order(colnames(comms_p.ca))]

# #combine
# allcomms <- bind_rows(comms.ca,comms_p.ca)
# allcomms <- allcomms[,c(1:32)] #only OG 25 right now
# 
# #dissimilarity
# bcdist.ca <- vegdist.ca(as.matrix(allcomms), method="bray")
# bcdist.ca <-as.matrix(bcdist.ca)

##2021
# within each year subset the data
div.ca21 <- comp.ca %>% filter(year=="2021")
div.ca21 <- div.ca21 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
div.ca21 <- div.ca21 %>% column_to_rownames("trt.b.y")
div.ca21 <- div.ca21[,c(16:54)]
# make comms.ca the long to add 4-letter codes
codes <- data.frame(sixletter = colnames(div.ca21))
codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
colnames(div.ca21) <- codes$sppcodes
colnames(div.ca21)[colnames(div.ca21) == "CACI"] <- "CAME" #change name of CACI to CAME
div.ca21 <- div.ca21[ order(row.names(div.ca21)), ]
# attach row of targets
comms_p.ca <- comms_p.ca[-c(99,205),] #remove "fd 50", not in comps data
all21 <- bind_rows(comms_p.ca, div.ca21)
all21 <- all21[,c(1:32)] #only OG 25 right now
#run div.ca or vegdiv.ca
div.camat <- diversity(as.matrix(all21))
div.camat21 <-data.frame(div.camat)
div.camat21$plot <- rownames(div.camat21)

##2022
# within each year subset the data
div.ca22 <- comp.ca %>% filter(year=="2022")
div.ca22 <- div.ca22 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
div.ca22 <- div.ca22 %>% column_to_rownames("trt.b.y")
div.ca22 <- div.ca22[,c(16:54)]
# make comms.ca the long to add 4-letter codes
#codes <- data.frame(sixletter = colnames(div.ca22))
#codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
colnames(div.ca22) <- codes$sppcodes
colnames(div.ca22)[colnames(div.ca22) == "CACI"] <- "CAME" #change name of CACI to CAME
div.ca22 <- div.ca22[ order(row.names(div.ca22)), ]
# attach row of targets
#comms_p.ca <- comms_p.ca[-c(99,205),] #remove "ir 50" and "fd 50", not in comps data
all22 <- bind_rows(comms_p.ca, div.ca22)
all22 <- all22[,c(1:32)] #only OG 25 right now
#run div.ca or vegdiv.ca
div.camat <- diversity(as.matrix(all22))
div.camat22 <-data.frame(div.camat)
div.camat22$plot <- rownames(div.camat22)

##2023
# within each year subset the data
div.ca23 <- comp.ca %>% filter(year=="2023")
div.ca23 <- div.ca23 %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=T)
div.ca23 <- div.ca23 %>% column_to_rownames("trt.b.y")
div.ca23 <- div.ca23[,c(16:54)]
# make comms.ca the long to add 4-letter codes
#codes <- data.frame(sixletter = colnames(div.ca23))
#codes <- codes %>% mutate(sppcodes = paste0(substr(sixletter, 1, 2), substr(sixletter, 4, 5))) # make 4 letter code to get cwm's
colnames(div.ca23) <- codes$sppcodes
colnames(div.ca23)[colnames(div.ca23) == "CACI"] <- "CAME" #change name of CACI to CAME
div.ca23 <- div.ca23[ order(row.names(div.ca23)), ]
# attach row of targets
#comms_p.ca <- comms_p.ca[-c(99,205),] #remove "ir 50" and "fd 50", not in comps data
all23 <- bind_rows(comms_p.ca, div.ca23)
all23 <- all23[,c(1:32)] #only OG 25 right now
#run div.ca or vegdiv.ca
div.camat <- diversity(as.matrix(all23))
div.camat23 <-data.frame(div.camat)
div.camat23$plot <- rownames(div.camat23)


## combine dissimilarity dataframe
div.ca <- bind_rows(div.camat21,div.camat22,div.camat23)
div.ca <- div.ca[,-c(1:208)]
div.ca <- separate(div.ca, plot, into = c("trt", "block", "subplot", "year"), sep = "\\.")
ggplot(div.ca, aes(x=div.camat,col=trt))+
  geom_boxplot()

ggplot(div.ca, aes(x=div.camat,y=full,col=trt))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~year, scales="free")

#save
write.csv(bcdiv.ca, "data/div_ca.csv" )


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
# bcdiv.ca <- vegdiv.ca(as.matrix(allcomms), method="bray")
# bcdiv.ca <-as.matrix(bcdiv.ca)

##2021
# within each year subset the data
div.wy21 <- comp.wy.wide %>% filter(year=="2021")
div.wy21 <- div.wy21 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
div.wy21 <- div.wy21 %>% column_to_rownames("trt.b.sub.y")
div.wy21 <- div.wy21[,c(7:62)]
# make comms.wy the long to add 4-letter codes
div.wy21 <- div.wy21[ order(row.names(div.wy21)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all21 <- bind_rows(comms_p.wy, div.wy21)
all21 <- all21[,c(1:25)] #only OG 25 right now
all21[is.na(all21)] <- 0
#run div.wy or vegdiv.wy
bcdiv.wymat <- vegdiv(as.matrix(all21),method = "bray")
bcdiv.wymat <-as.matrix(bcdiv.wymat)[c(1:256),]
bcdiv.wymat <-as.matrix(bcdiv.wymat)[,-c(1:256)]
#north first
filtered_dissimilarity_df <- data.frame(bcdiv.wymat) %>%
  select(matches("*\\.n\\.*"))
# save only pairwise between target 
bcdiv.wy21.n<- data.frame(
  div.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#south next
filtered_dissimilarity_df <- data.frame(bcdiv.wymat) %>%
  select(matches("*\\.s\\.*"))
# save only pairwise between target 
bcdiv.wy21.s<- data.frame(
  div.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#combine
bcdiv.wy21 <- bind_rows(bcdiv.wy21.n,bcdiv.wy21.s)

##2022
# within each year subset the data
div.wy22 <- comp.wy.wide %>% filter(year=="2022")
div.wy22 <- div.wy22 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
div.wy22 <- div.wy22 %>% column_to_rownames("trt.b.sub.y")
div.wy22 <- div.wy22[,c(7:62)]
# make comms.wy the long to add 4-letter codes
div.wy22 <- div.wy22[ order(row.names(div.wy22)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all22 <- bind_rows(comms_p.wy, div.wy22)
all22 <- all22[,c(1:25)] #only OG 25 right now
all22[is.na(all22)] <- 0
#run div.wy or vegdiv.wy
bcdiv.wymat <- vegdiv(as.matrix(all22),method = "bray")
bcdiv.wymat <-as.matrix(bcdiv.wymat)[c(1:256),]
bcdiv.wymat <-as.matrix(bcdiv.wymat)[,-c(1:256)]
#north first
filtered_dissimilarity_df <- data.frame(bcdiv.wymat) %>%
  select(matches("*\\.n\\.*"))
# save only pairwise between target 
bcdiv.wy22.n<- data.frame(
  div.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#south next
filtered_dissimilarity_df <- data.frame(bcdiv.wymat) %>%
  select(matches("*\\.s\\.*"))
# save only pairwise between target 
bcdiv.wy22.s<- data.frame(
  div.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#combine
bcdiv.wy22 <- bind_rows(bcdiv.wy22.n,bcdiv.wy22.s)

##2023
# within each year subset the data
div.wy23 <- comp.wy.wide %>% filter(year=="2023")
div.wy23 <- div.wy23 %>% unite(trt.b.sub.y, c(trt, block, subplot,year), sep = ".", remove=T)
div.wy23 <- div.wy23 %>% column_to_rownames("trt.b.sub.y")
div.wy23 <- div.wy23[,c(7:62)]
# make comms.wy the long to add 4-letter codes
div.wy23 <- div.wy23[ order(row.names(div.wy23)), ]
# attach row of targets
#comms_p.wy <- comms_p.wy[-c(99,205),] #remove "fd 50", not in comps data
all23 <- bind_rows(comms_p.wy, div.wy23)
all23 <- all23[,c(1:25)] #only OG 25 right now
all23[is.na(all23)] <- 0
#run div.wy or vegdiv.wy
bcdiv.wymat <- vegdiv(as.matrix(all23),method = "bray")
bcdiv.wymat <-as.matrix(bcdiv.wymat)[c(1:256),]
bcdiv.wymat <-as.matrix(bcdiv.wymat)[,-c(1:256)]
#north first
filtered_dissimilarity_df <- data.frame(bcdiv.wymat) %>%
  select(matches("*\\.n\\.*"))
# save only pairwise between target 
bcdiv.wy23.n<- data.frame(
  div.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#south next
filtered_dissimilarity_df <- data.frame(bcdiv.wymat) %>%
  select(matches("*\\.s\\.*"))
# save only pairwise between target 
bcdiv.wy23.s<- data.frame(
  div.wy=diag(as.matrix(filtered_dissimilarity_df)),
  id=colnames(filtered_dissimilarity_df))
#combine
bcdiv.wy23 <- bind_rows(bcdiv.wy23.n,bcdiv.wy23.s)

## combine dissimilarity dataframe
bcdiv.wy <- bind_rows(bcdiv.wy21,bcdiv.wy22,bcdiv.wy23)

#save
write.csv(bcdiv.wy, "data/bc_dissimilarity_wy.csv" )

#### figures 
#CA
x <- read.csv("data/comp_ca.csv")
x$trt <- tolower(x$trt)
x <- x %>% mutate(trt = str_replace(trt, "^r$", "rand")) #make r match rand in cwm df
droughtca <- x %>% unite(trt.b.y, c(trt, block, year), sep = ".", remove=F)
droughtca <- droughtca %>% select(trt.b.y,water)
bcdiv.ca <- merge(bcdiv.ca,droughtca, by.y="trt.b.y",by.x="id", all.x=T)
bcdiv.ca <- separate(bcdiv.ca, id, into = c("trt", "block", "year"), sep = "\\.")
bcdiv.ca$year <- as.factor(bcdiv.ca$year)
bcdiv.ca$block <- as.factor(bcdiv.ca$block)
bcdiv.ca$trt <- as.factor(bcdiv.ca$trt)
bcdiv.ca$water <- as.factor(bcdiv.ca$water)
bcdiv.ca <- bcdiv.ca %>% mutate(drought=as.factor(ifelse(water=="0.5","drt","cntl")))
droughtcolsca <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

summary(aov(div.ca~trt*drought*year, bcdiv.ca))
ggplot(bcdiv.ca, aes(y=div.ca, x=trt, fill=drought))+
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
bcdiv.wy <- merge(bcdiv.wy,droughtwy, by.y="trt.b.sub.y",by.x="id", all.x=T)
bcdiv.wy <- separate(bcdiv.wy, id, into = c("trt", "block", "subplot","year"), sep = "\\.")
bcdiv.wy$year <- as.factor(bcdiv.wy$year)
bcdiv.wy$block <- as.factor(bcdiv.wy$block)
bcdiv.wy$trt <- as.factor(bcdiv.wy$trt)
bcdiv.wy <- bcdiv.wy %>% mutate(drought=as.factor(ifelse(drought=="0","cntl","drt")))
#bcdiv.wy$drought <- as.factor(bcdiv.wy$drought)
droughtcolswy <- c("cntl"="skyblue", "drt"="tomato1") #create variable for color

summary(aov(div.wy~trt*drought*year, bcdiv.wy))
ggplot(bcdiv.wy, aes(y=div.wy, x=trt, fill=drought))+
  geom_boxplot()+
  #geom_smooth(method="lm")+
  scale_fill_manual(values=droughtcolswy)+
  facet_wrap(~year, scales="fixed")+
  labs(x=" ",y="Dissimilarity (bray-curtis)")+ #, fill="drought treatment")+
  theme_classic()