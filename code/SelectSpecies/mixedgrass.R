########### Mixed-grass restoration experimental assemblages
########### Daniel Laughlin

setwd("C:/Users/jesse/Desktop/GitHub/happygrass/code/daniels_code")

library(FD)
library(mice)
library(RColorBrewer)
library(rlist)
library(labdsv)
library(vegan)
library(sads)
library(erer)

# must load newSelectSpecies function
source("newSelectSpecies.R")

traits <- read.csv("mixedgrass.csv", header=TRUE, row.names=1)
traits <- traits[traits$use==1,] # subset use=1
traits$PLSg.m2.mono <- traits$PLSlbs.acre.mono * (453.59237 / 4046.8564224) #convert lb/acre to g/m2
#traits$PLSg.plot.mono <- traits$PLSg.m2.mono * (1.275/1000) #convert g/m2 to g/plot
#traits$seedmass.mg <- 453592.37 * (1/traits$seeds.lbs)
#traits$seeds.m2 <- traits$PLSmg.m2.mono / traits$seedmass.mg

#scale traits
traits$srl = scale(log(traits$srl))
traits$ldmc = scale(log(traits$ldmc))
traits$leafn = scale(log(traits$leafn))
traits$lop = scale(traits$lop)
traits$rootdiam = scale(log(traits$rootdiam))
traits$sla = scale(log(traits$sla))
traits$rdmc = scale(log(traits$rdmc))

### pca
pca <- princomp(traits[,3:9], cor=TRUE)
summary(pca)
biplot(pca)
traits$pc1 <- pca$scores[,1]
traits$pc2 <- pca$scores[,2]

#set color scheme
colors = brewer.pal(12,"Paired")
traits$cols = c(colorRampPalette(colors)(nrow(traits)))

### Subset monocots and dicots
grams <- subset(traits, graminoid==1)
forbs <- subset(traits, graminoid==0)


set.seed(20021924) # EC Pielou's birthday
N=64
###### Create invasion resistant assemblages ######
ir.assemblages <- list()
ir.output <- list()
for(i in 1:N){
  ### Randomly sample 6 monocots and 9 dicots for each random species pool
  sub.grams <- grams[sample(nrow(grams),6,replace=FALSE),]
  sub.forbs <- forbs[sample(nrow(forbs),9,replace=FALSE),]
  pool <- rbind(sub.grams, sub.forbs)
  
  #### trait matrices
  traits.ir <- as.matrix(pool[,c("leafn","srl")])
  colnames(traits.ir) <- c("leafn","srl")
  rownames(traits.ir) <- rownames(traits.ir)
  traits.fd <- as.matrix(pool[,c("veg")])
  rownames(traits.fd) <- rownames(traits.ir)
  colnames(traits.fd) <- c("veg")
  
  ######## trait constraints # low Leaf N, high graminoid cover
  leafn.low = quantile(traits$leafn, 0.25)
  srl.high = quantile(traits$srl, 0.7557) #many digits to reproduce innocuous bug "quantile(traits$leafn, 0.75)"
  constraints.ir <- as.matrix(t(c(leafn.low, srl.high)))
  colnames(constraints.ir) <- c("leafn.low","srl.high")
  rownames(constraints.ir) <- c("constraints")
  
  ### Invasion resistant assemblages
  m1 <- newSelectSpecies(t2c=traits.ir, constraints=constraints.ir, t2d=traits.fd, obj="QH", capd=FALSE)
  if(m1$convergence == 0) {
    m1.out <- data.frame(trt=paste("ir",i), species=rownames(pool), prob=round(m1$prob,6), cols=pool$cols,
                         graminoid=pool$graminoid, c4=pool$c4, veg=pool$veg, g.m2=pool$PLSg.m2.mono,
                         seeds.g=pool$seeds.lbs/453.592, pls.mult=pool$pls.multiplier,
                         pc1=pool$pc1, pc2=pool$pc2)
    m1.out$seedrate.g.plot <- 3* round(m1.out$prob * (m1.out$g.m2 * (1.275) ) * m1.out$pls.mult, 4) # convert g/m2.mono to g/plot(1.275m2)
    m1.sort <- m1.out[order(-m1.out$prob),]
    ir.assemblages[[i]] <- m1.sort
  }else{
    ir.assemblages[[i]] <- NA
  }
  ir.output[[i]] <- m1
}
conv.ir <- (100 - sum(is.na(ir.assemblages))) / 100 # proportion converged


###### Create drought tolerant assemblages ######
dt.assemblages <- list()
dt.output <- list()
for(i in 1:N){
  
  ### Randomly sample 6 monocots and 9 dicots for each random species pool
  sub.grams <- grams[sample(nrow(grams),6,replace=FALSE),]
  sub.forbs <- forbs[sample(nrow(forbs),9,replace=FALSE),]
  pool <- rbind(sub.grams, sub.forbs)
  
  #### trait matrix
  traits.dt <- as.matrix(pool[,c("lop", "ldmc")])
  colnames(traits.dt) <- c("lop", "ldmc")
  rownames(traits.dt) <- rownames(traits.dt)
  traits.fd <- as.matrix(as.numeric(pool[,c("rootdiam")]))
  rownames(traits.fd) <- rownames(traits.dt)
  colnames(traits.fd) <- "rootdiam"
  
  ######## trait constraints high LDMC, low LOP, high c4
  ldmc.high = quantile(traits$ldmc, 0.75)
  lop.low = quantile(traits$lop, 0.25)
  constraints.dt <- as.matrix(t(c(lop.low, ldmc.high))) ### constraints for drought tolerance
  colnames(constraints.dt) <- c("lop", "ldmc")
  
  ### Drought tolerance
  m2 <- newSelectSpecies(t2c=traits.dt, constraints=constraints.dt, t2d=traits.fd, obj="QH", capd=FALSE)
  if(m2$convergence == 0) {
    m2.out <- data.frame(trt=paste("dt",i), species=rownames(pool), prob=round(m2$prob,6), cols=pool$cols,
                         graminoid=pool$graminoid, c4=pool$c4, veg=pool$veg, g.m2=pool$PLSg.m2.mono,
                         seeds.g=pool$seeds.lbs/453.592, pls.mult=pool$pls.multiplier,
                         pc1=pool$pc1, pc2=pool$pc2)
    m2.out$seedrate.g.plot <- 3* round(m2.out$prob * (m2.out$g.m2 * (1.275) ) * m2.out$pls.mult, 4) # convert g/m2 to g/plot(1.275m2)
    m2.sort <- m2.out[order(-m2.out$prob),]
    dt.assemblages[[i]] <- m2.sort
  }else{
    dt.assemblages[[i]] <- NA
  }
  dt.output[[i]] <- m2
}
conv.dt <- (100 - sum(is.na(dt.assemblages))) / 100 # proportion converged


###### Create functionally diverse assemblages ######
fd.assemblages <- list()
fd.output <- list()
for(i in 1:N){
  
  ### Randomly sample 6 monocots and 9 dicots for each random species pool
  sub.grams <- grams[sample(nrow(grams),6,replace=FALSE),]
  sub.forbs <- forbs[sample(nrow(forbs),9,replace=FALSE),]
  pool <- rbind(sub.grams, sub.forbs)
  
  #### trait matrix
  traits.fd <- as.matrix(pool[,c("leafn","srl","ldmc","lop","rootdiam","veg")])
  
  ### Functional diversity
  m3 <- newSelectSpecies(t2c=traits.fd, constraints=NULL, t2d=traits.fd, obj="QH", capd=FALSE)
  if(m2$convergence == 0) {
    m3.out <- data.frame(trt=paste("fd",i), species=rownames(pool), prob=round(m3$prob,6), cols=pool$cols,
                         graminoid=pool$graminoid, c4=pool$c4, veg=pool$veg, g.m2=pool$PLSg.m2.mono,
                         seeds.g=pool$seeds.lbs/453.592, pls.mult=pool$pls.multiplier,
                         pc1=pool$pc1, pc2=pool$pc2)
    m3.out$seedrate.g.plot <- 3* round(m3.out$prob * (m3.out$g.m2 * (1.275) ) * m3.out$pls.mult, 4) # convert g/m2 to g/plot(1.275m2)
    m3.sort <- m3.out[order(-m3.out$prob),]
    fd.assemblages[[i]] <- m3.sort
  }else{
    fd.assemblages[[i]] <- NA
  }
  fd.output[[i]] <- m3
}
conv.fd <- (100 - sum(is.na(fd.assemblages))) / 100 # proportion converged


###### Create random assemblages #########
### cite Ulrich et al. 2010 Oikos and Harpole and Tilman 2006 Ecol Lett for lognormal in grasslands
random.assemblages <- list()
for (i in 1:N){
  sub.grams <- grams[sample(nrow(grams),6,replace=FALSE),]
  sub.forbs <- forbs[sample(nrow(forbs),9,replace=FALSE),]
  pool <- rbind(sub.grams, sub.forbs)
  sim.abun <- sads::rsad(S = 15, frac=0.1, sad="lnorm", coef=list(meanlog=4.2, sdlog=1.8), zeroes=TRUE)
  random.abun <- sim.abun/sum(sim.abun)
  assemblage <- data.frame(trt=paste("r",i), species = rownames(pool), prob = random.abun, cols=pool$cols, 
                           graminoid=pool$graminoid, c4=pool$c4, veg=pool$veg, g.m2=pool$PLSg.m2.mono,
                           seeds.g=pool$seeds.lbs/453.592, pls.mult=pool$pls.multiplier,
                           pc1=pool$pc1, pc2=pool$pc2)
  assemblage$seedrate.g.plot <- 3* round(assemblage$prob * (assemblage$g.m2 * (1.275) ) * assemblage$pls.mult, 4) # convert g/m2 to g/plot(1.275m2)
  rownames(assemblage) <- rownames(pool)
  assemblage.sorted <- assemblage[order(-assemblage$prob),]
  random.assemblages[[i]] <- assemblage.sorted
}
conv.rand <- (100 - sum(is.na(random.assemblages))) / 100 # proportion converged

### Test for incompatible trait profile
try(if(conv.ir < 1) stop("incompatible 'ir' trait profile"))
try(if(conv.dt < 1) stop("incompatible 'dt' trait profile"))
try(if(conv.fd < 1) stop("incompatible 'fd' trait profile"))
try(if(conv.rand < 1) stop("incompatible 'random' trait profile"))

#### End of program



### plot NMDS of assemblages
ir.df <- data.frame(list.stack(ir.assemblages))
dt.df <- data.frame(list.stack(dt.assemblages))
fd.df <- data.frame(list.stack(fd.assemblages))
rand.df <- data.frame(list.stack(random.assemblages))
comms <- rbind(ir.df, dt.df, fd.df, rand.df)
comms <- rbind(fd.df, rand.df)
comms <- matrify(comms[,1:3])
comms <- comms[,order(colnames(comms))]
groups <- c(rep("ir",N), rep("dt",N), rep("fd",N), rep("rand",N))
groups <- c(rep("fd",N), rep("rand",N))
nms <- metaMDS(comms)
plot(nms, type="t")
ordiellipse(nms, groups, col=1:4, conf=0.95)
trait.matrix <- traits[order(rownames(traits)),]
trait.matrix <- as.matrix(trait.matrix[,2:11])
comm.matrix <- as.matrix(comms)
cwm <- FD::functcomp(as.matrix(trait.matrix), as.matrix(comms), bin.num=c("graminoid","veg","c4"))
nms2 <- metaMDS(cwm)
plot(nms2, type="t")
ordiellipse(nms2, groups, col=1:4, conf=0.95)
vectors <- envfit(nms2,cwm)
plot(vectors,p.max=0.05)




### Test for graminoid proportion
gram.prop <- c()
for(i in 1:N) gram.prop[i] <- ir.assemblages[[i]]$prob %*% ir.assemblages[[i]]$graminoid
ir.gram <- c(mean=mean(gram.prop), range=range(gram.prop))
for(i in 1:N) gram.prop[i] <- dt.assemblages[[i]]$prob %*% dt.assemblages[[i]]$graminoid
dt.gram <- c(mean=mean(gram.prop), range=range(gram.prop))
for(i in 1:N) gram.prop[i] <- fd.assemblages[[i]]$prob %*% fd.assemblages[[i]]$graminoid
fd.gram <- c(mean=mean(gram.prop), range=range(gram.prop))
for(i in 1:N) gram.prop[i] <- random.assemblages[[i]]$prob %*% random.assemblages[[i]]$graminoid
random.gram <- c(mean=mean(gram.prop), range=range(gram.prop))
gram.df <- rbind(ir.gram, dt.gram, fd.gram, random.gram)

### Test for c4 proportion
c4.prop <- c()
for(i in 1:N) c4.prop[i] <- ir.assemblages[[i]]$prob %*% ir.assemblages[[i]]$c4
ir.c4 <- c(mean=mean(c4.prop), range=range(c4.prop))
for(i in 1:N) c4.prop[i] <- dt.assemblages[[i]]$prob %*% dt.assemblages[[i]]$c4
dt.c4 <- c(mean=mean(c4.prop), range=range(c4.prop))
for(i in 1:N) c4.prop[i] <- fd.assemblages[[i]]$prob %*% fd.assemblages[[i]]$c4
fd.c4 <- c(mean=mean(c4.prop), range=range(c4.prop))
for(i in 1:N) c4.prop[i] <- random.assemblages[[i]]$prob %*% random.assemblages[[i]]$c4
random.c4 <- c(mean=mean(c4.prop), range=range(c4.prop))
c4.df <- rbind(ir.c4, dt.c4, fd.c4, random.c4)

### Test for veg proportion
veg.prop <- c()
for(i in 1:N) veg.prop[i] <- ir.assemblages[[i]]$prob %*% ir.assemblages[[i]]$veg
ir.veg <- c(mean=mean(veg.prop), range=range(veg.prop))
for(i in 1:N) veg.prop[i] <- dt.assemblages[[i]]$prob %*% dt.assemblages[[i]]$veg
dt.veg <- c(mean=mean(veg.prop), range=range(veg.prop))
for(i in 1:N) veg.prop[i] <- fd.assemblages[[i]]$prob %*% fd.assemblages[[i]]$veg
fd.veg <- c(mean=mean(veg.prop), range=range(veg.prop))
for(i in 1:N) veg.prop[i] <- random.assemblages[[i]]$prob %*% random.assemblages[[i]]$veg
random.veg <- c(mean=mean(veg.prop), range=range(veg.prop))
veg.df <- rbind(ir.veg, dt.veg, fd.veg, random.veg)

gram.df
c4.df
veg.df



##### Save output into txt files
write.list(ir.assemblages, file = "ir.assemblages.txt", row.names = TRUE)
ir.assemblages.csv <- read.csv("ir.assemblages.txt", sep=",")
ir.assemblages.print <- ir.assemblages.csv[,c(3,4,5,7,13)]
write.table(ir.assemblages.print, "ir.assemblages.print.csv", sep=",")
#capture.output(ir.output, file = "ir.output.txt")

write.list(dt.assemblages, file = "dt.assemblages.txt", row.names = TRUE)
dt.assemblages.csv <- read.csv("dt.assemblages.txt", sep=",")
dt.assemblages.print <- dt.assemblages.csv[,c(3,4,5,7,13)]
write.table(dt.assemblages.print, "dt.assemblages.print", sep=",")
#capture.output(dt.output, file = "dt.output.txt")

write.list(fd.assemblages, file = "fd.assemblages.txt", row.names = TRUE)
fd.assemblages.csv <- read.csv("fd.assemblages.txt", sep=",")
fd.assemblages.print <- fd.assemblages.csv[,c(3,4,5,7,13)]
write.table(fd.assemblages.print, "fd.assemblages.print", sep=",")
#capture.output(fd.output, file = "fd.output.txt")

write.list(random.assemblages, file = "random.assemblages.txt", row.names = TRUE)
random.assemblages.csv <- read.csv("random.assemblages.txt", sep=",")
random.assemblages.print <- random.assemblages.csv[,c(3,4,5,7,13)]
write.table(random.assemblages.print, "random.assemblages.print", sep=",")


###### check if seed rates are good and seed supply #######
### Aim for > 1500 pure live seed PLS / m2 #Barr et al. 2017 Rest Ecol
### Aim for > 15 kg/ha (1.5 g/m2)  https://www.npss.sk.ca/docs/2_pdf/Seeding_Rate_Conversion.pdf
### Aim for > 12 lbs PLS/acre = 1.34 g/m2 #https://webapp.agron.ksu.edu/agr_social/m_eu_article.throck?article_id=1157

allplots <- rbind(ir.assemblages.csv, dt.assemblages.csv, fd.assemblages.csv, random.assemblages.csv)
allplots <- allplots[!allplots$species=="species",]
allplots <- allplots[!allplots$species=="",]
sp.seed <- allplots[,c(4,15)]
sp.seed[,2] <- as.numeric(as.character(sp.seed[,2]))
seed.needed <- aggregate(data=sp.seed, seedrate.g.plot ~ species, FUN="sum")
colnames(seed.needed) <- c("species","total.seed.needed")
total.seed <- traits[,c("spcode4","totalseed.g")]
colnames(total.seed) <- c("species","totalseed.g")

### check seed supply # should be FALSE
seed.total.comp <- merge(seed.needed, total.seed, by="species")
seed.total.comp$factor <- seed.total.comp$totalseed.g / seed.total.comp$total.seed.needed
try(any(seed.total.comp$factor < 1))

### check number of seeds sown per plot # mean should be > 1500 seed/m2
allplots$seed.count <- as.numeric(as.character(allplots$seeds.g)) *
                                    as.numeric(as.character(allplots$seedrate.g.plot))
seed.per.plot <- aggregate(data=allplots, seed.count ~ trt, FUN="sum")
mean(seed.per.plot[,2]); range(seed.per.plot[,2]); hist(seed.per.plot[,2])

### check seed mass sown per plot # mean should be > 3 g/m2
plot.seed <- allplots[,c(3,15)]
plot.seed[,2] <- as.numeric(as.character(plot.seed[,2]))
total.seed.rate <- aggregate(data=plot.seed, seedrate.g.plot ~ trt, FUN="sum")
mean(total.seed.rate[,2]); range(total.seed.rate[,2]); hist(total.seed.rate[,2])

######   Final print out   ######
allplots <- tidyr::separate(allplots, col=Result, into=c("result","block"))
allplots <- allplots[order(as.numeric(allplots$block)),]
allplots <- allplots[,c("block","trt","species","prob","graminoid","seedrate.g.plot")]
#write.table(allplots, "allplot.assemblages.csv", sep=",", row.names=FALSE)




###### Winter seeding rate ######
fall <- allplots
fall$seedrate.g.plot <- as.numeric(as.character(fall$seedrate.g.plot)) #make numeric
fall$prob <- as.numeric(as.character(fall$prob)) #make numeric
fall$cov[1:15] <- c(0.15,0.05,0.25,0.02,0.05,0.02,0,0,0,0.02,0.01,0.05,0.05,0,0) #dummy data
fall$diffratio <- fall$cov/fall$prob
fall$fallrate.g.plot <- fall$seedrate.g.plot
for(i in 1:nrow(fall)){
  ifelse(fall$diffratio[i] > 1, fall$fallrate.g.plot[i] <- 0, 
         fall$fallrate.g.plot[i] <- ((1-fall$diff[i])*fall$seedrate.g.plot[i]))
}

#fall$fallrate.g.plot[fall$fallrate.g.plot < 0.15 ] <- 0.15
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "ACMI"] <- 0.1
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "ARFR"] <- 0.1
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "ARLU"] <- 0.1
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "KOMA"] <- 0.1
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "POSE"] <- 0.1
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "POPE"] <- 0.1
#fall$fallrate.g.plot[fall$fallrate.g.plot == 0.1500 & fall$species == "SPCR"] <- 0.1

#total.fall.seed.rate <- aggregate(data=fall, fallrate.g.plot ~ trt, FUN="sum")
#mean(total.fall.seed.rate[,2]); range(total.fall.seed.rate[,2]); hist(total.fall.seed.rate[,2])





#####   Growing plugs for outplanting   ####
plugs <- allplots
plugs$cumsum <- 0
plugs$prob <- as.numeric(as.character(plugs$prob)) #make numeric
for (i in levels(plugs$trt)){
  plugs$cumsum[plugs$trt == i] <- cumsum(plugs$prob[plugs$trt == i])
} # compute cumsum
plugs$seedlings <- ceiling(round(plugs$prob,1)*14) #round to one decimal and make 16 plants, then round up to integer
plugs$seedlings[plugs$cumsum > 0.83] <- 0 #no seedlings if spp does not contribute to 80% biomass
plugs$seedlings[plugs$seedlings == 0 & plugs$cumsum < 0.83] <- 1 #if spp contributes to 80% biomass and >0.05 prob then add 1 plants
plugs$seedlings[plugs$prob == plugs$cumsum & plugs$seedlings == 0] <- ceiling(round(plugs$prob[plugs$prob == plugs$cumsum & plugs$seedlings == 0],1)*14)
plugs$seedlings[plugs$seedlings == 0 & plugs$prob > 0.1] <- 2
plugs[plugs$prob > 0.1 & plugs$seedlings == 0,] # check if common species were inadvertantly set to zero
sum(plugs$seedlings)
plugsums <- aggregate(plugs$seedlings, by=list(species=plugs$species), FUN=sum)
write.table(plugs, "plugs.csv", sep=",", row.names=FALSE)
write.table(plugsums, "plugsums.csv", sep=",", row.names=FALSE)


##### show PCA biplots with circles proportional to species probs
for (i in 1:length(ir.assemblages)){
plot(traits$pc1, traits$pc2, col="white", xlab="PC1", ylab="PC2")
text(traits$pc1, traits$pc2, rownames(traits) )
# loadings
x0 <- c(0, 0, 0, 0)
y0 <- c(0, 0, 0, 0)
x1 <- pca$loadings[,1]*3
y1 <- pca$loadings[,2]*3
arrows(x0, y0, x1, y1, col = 4, lwd = 2, d = 3,  ticktype = "detailed", add=TRUE)
text(x1*1.1, y1*1.1, colnames(traits)[3:9], col=4, add=TRUE)
points(ir.assemblages[[i]]$pc1, ir.assemblages[[i]]$pc2, cex=ir.assemblages[[i]]$prob*100)
}

for (i in 1:length(dt.assemblages)){
plot(traits$pc1, traits$pc2, col="white", xlab="PC1", ylab="PC2")
text(traits$pc1, traits$pc2, rownames(traits) )
# loadings
x0 <- c(0, 0, 0, 0)
y0 <- c(0, 0, 0, 0)
x1 <- pca$loadings[,1]*3
y1 <- pca$loadings[,2]*3
arrows(x0, y0, x1, y1, col = 4, lwd = 2, d = 3,  ticktype = "detailed", add=TRUE)
text(x1*1.1, y1*1.1, colnames(traits)[3:9], col=4, add=TRUE)
points(dt.assemblages[[i]]$pc1, dt.assemblages[[i]]$pc2, cex=dt.assemblages[[i]]$prob*100)
}

for (i in 1:length(fd.assemblages)){
plot(traits$pc1, traits$pc2, col="white", xlab="PC1", ylab="PC2")
text(traits$pc1, traits$pc2, rownames(traits) )
# loadings
x0 <- c(0, 0, 0, 0)
y0 <- c(0, 0, 0, 0)
x1 <- pca$loadings[,1]*3
y1 <- pca$loadings[,2]*3
arrows(x0, y0, x1, y1, col = 4, lwd = 2, d = 3,  ticktype = "detailed", add=TRUE)
text(x1*1.1, y1*1.1, colnames(traits)[3:9], col=4, add=TRUE)
points(fd.assemblages[[i]]$pc1, fd.assemblages[[i]]$pc2, cex=fd.assemblages[[i]]$prob*100)
}

for (i in 1:length(random.assemblages)){
plot(traits$pc1, traits$pc2, col="white", xlab="PC1", ylab="PC2")
text(traits$pc1, traits$pc2, rownames(traits) )
# loadings
x0 <- c(0, 0, 0, 0)
y0 <- c(0, 0, 0, 0)
x1 <- pca$loadings[,1]*3
y1 <- pca$loadings[,2]*3
arrows(x0, y0, x1, y1, col = 4, lwd = 2, d = 3,  ticktype = "detailed", add=TRUE)
text(x1*1.1, y1*1.1, colnames(traits)[3:9], col=4, add=TRUE)
points(random.assemblages[[i]]$pc1, random.assemblages[[i]]$pc2, cex=random.assemblages[[i]]$prob*100)
}







### Test SADs in mixed grass prairie
# prepare data
hpg <- read.csv("HPGRS_basal_foliar_cov_2003to2015.csv", header=TRUE)
hpgmeta <- hpg[,1:8]
hpgcomm <- hpg[,grepl("_abs_foliar",names(hpg))]
hpgcomm <- hpgcomm[, order(colnames(hpgcomm))]
rownames(hpgcomm) <- seq(1:nrow(hpgcomm))
hpg.long <- labdsv::dematrify(hpgcomm, thresh=0)

# fit SADs
m <- list()
meanlog <- c()
sdlog <- c()
for(i in 1:nrow(hpgcomm)){
  m[[i]] <- fitsad(x=hpg.long[hpg.long$sample==i,]$abundance*100, sad="lnorm")
  meanlog[i] <- coef(m[[i]])[1]
  sdlog[i] <- coef(m[[i]])[2]
}

#mean parameters across
hist(meanlog);mean(meanlog)
hist(sdlog);mean(sdlog)




### Figure explaining experimental design

#alternative color scheme to signify grass and forbs better (not working yet)
#colors = brewer.pal(6:11,"Spectral")
#plot(NULL, xlim=c(0,length(colors)), ylim=c(0,1), 
#     xlab="", ylab="", xaxt="n", yaxt="n")
#rect(0:(length(colors)-1), 0, 1:length(colors), 1, col=colors)
#traits$cols.new[traits$graminoid==0] = c(colorRampPalette(colors[1:6])(nrow(traits[traits$graminoid==0,])))
#traits$cols.new[traits$graminoid==1] = c(colorRampPalette(colors[7:11])(nrow(traits[traits$graminoid==1,])))

### Plot of species pool
par(mfrow=c(1,1))
plot(-1:1,-1:1, col="white",xlab="",ylab="", xaxt='n', yaxt='n', bty='n')
legend("topleft",legend=traits$species[traits$graminoid==1], fill=traits$cols[traits$graminoid==1],
       text.font=3, bty="n", title="Grasses", cex=2)
legend("topright",legend=traits$species[traits$graminoid==0], fill=traits$cols[traits$graminoid==0],
       text.font=3, bty="n", title="Forbs", cex=2)

par(mfrow=c(4,3),mar=c(0,0,0,0))
pie(dt.assemblages[[1]]$prob, labels=NA, col=as.character(dt.assemblages[[1]]$cols))
pie(dt.assemblages[[2]]$prob, labels=NA, col=as.character(dt.assemblages[[2]]$cols))
pie(dt.assemblages[[4]]$prob, labels=NA, col=as.character(dt.assemblages[[4]]$cols))

pie(ir.assemblages[[1]]$prob, labels=NA, col=as.character(ir.assemblages[[1]]$cols))
pie(ir.assemblages[[2]]$prob, labels=NA, col=as.character(ir.assemblages[[2]]$cols))
pie(ir.assemblages[[3]]$prob, labels=NA, col=as.character(ir.assemblages[[3]]$cols))

pie(fd.assemblages[[1]]$prob, labels=NA, col=as.character(fd.assemblages[[1]]$cols))
pie(fd.assemblages[[2]]$prob, labels=NA, col=as.character(fd.assemblages[[2]]$cols))
pie(fd.assemblages[[3]]$prob, labels=NA, col=as.character(fd.assemblages[[3]]$cols))

pie(random.assemblages[[4]]$prob, labels=NA, col=as.character(random.assemblages[[4]]$cols))
pie(random.assemblages[[6]]$prob, labels=NA, col=as.character(random.assemblages[[6]]$cols))
pie(random.assemblages[[7]]$prob, labels=NA, col=as.character(random.assemblages[[7]]$cols))



### phylogenetics
library(brranching)
traits$species <- gsub("Pascopyrum smithii","Elymus smithii",traits$species)
tree <- phylomatic(traits$species, get="POST", storedtree="zanne2014" )
tree$tip.label <- Hmisc::capitalize(gsub("_"," ",tree$tip.label))
tree$tip.label <- gsub("Elymus smithii","Pascopyrum smithii",tree$tip.label)
plot(tree)

###pca of species in trait space
traits$species <- gsub("Elymus smithii","Pascopyrum smithii",traits$species)
traits$spcode <- rownames(traits)
rownames(traits) <- traits$species
p.dist.mat <- cophenetic(tree)         ### trick to ensure correct order
traits <- traits[row.names(p.dist.mat),] ### trick to ensure correct order

library(phytools)
pca <- phyl.pca(tree, traits[,3:9], mode="cov")
summary(pca)
rownames(traits) <- traits$spcode
biplot(pca, xlabs=traits$spcode)
traits$pc1 <- pca$scores[,1]
traits$pc2 <- pca$scores[,2]

