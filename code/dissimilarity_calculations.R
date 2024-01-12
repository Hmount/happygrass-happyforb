
dat <- read.csv("data/cwm_wy.csv")
dat$year <- as.factor(dat$year)
dat$block <- as.factor(dat$block)
dat$trt <- as.factor(dat$trt)
dat$subplot <- as.factor(dat$subplot)
FDdat <- read.csv("data/cwm_raoq_wy.csv")
FDdat$year <- as.factor(FDdat$year)
FDdat$block <- as.factor(FDdat$block)
FDdat$trt <- as.factor(FDdat$trt)
FDdat$subplot <- as.factor(FDdat$subplot)
## Calculating dissimilarity/ distance (euclidian) from targets
#select relevent CWM traits per DT or IR, rows are communities
distpred <- dat %>% filter(year=="0")#filter preds
dat <- dat %>% filter(year!="0")#filter preds
FDdistpred <- FDdat %>% filter(year=="0")#filter preds
FDdat <- FDdat %>% filter(year!="0")#filter preds
dist21.wy <- dat %>% select(c(block,trt,subplot,year,leafn,srl))
dist21.wy <- merge(dist21.wy,FDdat[,c(2,4:6,12)], all.x = T)#, by=c(block,trt,year,subplot))# add veg FD data
dist21.wy2 <- distpred %>% select(c(block,trt,year,leafn,srl))
dist21.wy2 <- merge(dist21.wy2,FDdistpred[,c(2,4:6)], all = T)# add veg FD data
dist21.wy3 <- merge(dist21.wy,dist21.wy2, all.x = T)#, by=c("trt", "block", "year"))
#dist21.wy3 <- bind_rows(dist21.wy,dist21.wy2)#, all.x = T)#, by=c("trt", "block", "year"))
dist21.wy3 <- dist21.wy3 %>% 
  unite(trt.b.sub.y, c(trt, block, subplot, year), sep = ".", remove=T) # make unique plot variable
dist21.wy3 <- dist21.wy3 %>% column_to_rownames("trt.b.sub.y")

#mnake into matrix
#run dist or vegdist
testdistmat <- dist(as.matrix(dist21.wy3),method = "euclidean")#,diag=T)

#make save-able dataframe of pairwise crosses of interest

