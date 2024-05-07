#### zero inflated beta regression
library(tidyverse)
library(glmmTMB)

## CA
# data

# find relativized cover in invasive abundance 
suballca <- suballca %>% mutate(invcov.rel = FESPER/(native.cover+inv.grass.cov))
mod1 <- glmmTMB(invcov.rel ~ drought * distir + (1|structure) + (1|trt),
        ziformula= ~ trt * drought * distir + (1|structure),
        data = suballca, family=beta_family())
summary(mod1)

#fix non-convergence
#step one: determine the effect of coefficents on the zero inflation componant
order(fixef(mod1)$zi)
#examine diagnostic plots
plot(simulateResiduals(mod1))

ttt <- suballca %>% filter(!is.na(log.invg))
hist(ttt$log.invg)
ggplot(suballca, aes(x = distdt, y = invcov.rel, col=trt)) +
  geom_point()+
  geom_smooth(method="lm")

quantile(suballca$native.cover,.05) #find lower .05
quantile(suballca$distir,.95) #find upper .95
library(ggeffects)
x <- ggpredict(mod1,c("trt", "distir [c(.72,2.51)]","drought","native.cover [c(.17,1.42)]")) #all smooths lines
plot(x)
distmod <- lmer(invcov.rel ~ distir * drought + (1 | structure), data = suballca)
summary(distmod) 
anova(distmod) #only drought important (Rdiam was, but not really relevant)
emm.ca <- emmeans(m2rd.ca, c("rootdiam","drought"))
pairs(emm.ca)
#view
library(ggpubr)
dtca.23<-ggplot(suballca, aes(y=invcov.rel,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  stat_cor(label.y = c(.2,.22))+
  labs(y="relative cover invasive grass")+
  #stat_regline_equation(label.x = 2.5) +
  #  annotate("text", x = min(data$x), y = max(data$y), 
  #           label = paste("R^2 =", round(r_squared, 3), "p-value =", round(p_value, 3)),
  #           hjust = 0, vjust = 1)
  # #  stat_regline_equation(aes(label = paste(after_stat(eq.label), 
  #                      paste("R^2 =", format(summary(model)$r.squared, digits = 2)))))+#, 
  #                           # sep = "*\", \"p-value =", format(summary(model)$coefficients[8], digits = 3), sep = "*, \""))))+
  # facet_wrap(~year)+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
irca.23<-ggplot(suballca, aes(y=invcov.rel,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="relative cover invasive grass")+
  stat_cor(label.y = c(.24,.22))+
  theme_ggeffects()+
  theme(legend.position = "none")
fdca.23<-ggplot(suballca, aes(y=invcov.rel,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="relative cover invasive grass")+
  stat_cor()+
  theme_ggeffects()+
  theme(legend.position = "none")
distplotsca<-ggarrange(irca,dtca,fdca, nrow=2, ncol=2,common.legend = T )

#plots 2/8
boxca.23 <-ggplot(suballca, aes(y=invcov.rel ,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  #facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="relative cover invasive grass", fill="seed trt")+
  theme_ggeffects()

ggarrange(irca,dtca,fdca, boxca,nrow=2, ncol=2,common.legend = F)
annotate_figure(ggarrange(irca.23,dtca.23,fdca.23, boxca.23,nrow=2, ncol=2,common.legend = F), "2023")

## WY model
#prepare var
suballwy <- suballwy %>% mutate(brte.cov.rel = BRTE/totcov.plot) #relative proportion of BRTE (not proportion of all invasive (yet))
#model
mod2 <- glmmTMB(brte.cov.rel ~ trt * drought * distir + (1|block),
                ziformula= ~ trt * drought * distir + (1|block), #does not run with distir???
                data = suballwy, family=beta_family())
summary(mod2)
hist(suballwy$brte.cov.rel)

#fix non-convergence
#step one: determine the effect of coefficents on the zero inflation componant
order(fixef(mod2)$zi)
#examine diagnostic plots
plot(simulateResiduals(mod2))

ggplot(suballwy, aes(x = nativecov.plot, y = invcov.rel)) +
  geom_point()
ggplot(suballwy, aes(x = distir, y = BRTE, col=trt)) +
  geom_point()+
  geom_smooth(method="lm")

dtwy<-ggplot(suballwy, aes(y=BRTE,x=distdt,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  stat_cor(label.x = 1.85)+
  labs(y="relative cover BRTE")+
  #stat_regline_equation(label.x = 2.5) +
  #  annotate("text", x = min(data$x), y = max(data$y), 
  #           label = paste("R^2 =", round(r_squared, 3), "p-value =", round(p_value, 3)),
  #           hjust = 0, vjust = 1)
  # #  stat_regline_equation(aes(label = paste(after_stat(eq.label), 
  #                      paste("R^2 =", format(summary(model)$r.squared, digits = 2)))))+#, 
  #                           # sep = "*\", \"p-value =", format(summary(model)$coefficients[8], digits = 3), sep = "*, \""))))+
  # facet_wrap(~year)+
  #geom_hline(yintercept =0,col="black")+
  theme_ggeffects()
irwy<-ggplot(suballwy, aes(y=BRTE,x=distir,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="relative cover BRTE")+
  stat_cor()+
  theme_ggeffects()+
  theme(legend.position = "none")
fdwy<-ggplot(suballwy, aes(y=BRTE,x=distfd,color=drought))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_manual(values=droughtcols)+
  labs(y="relative cover BRTE")+
  stat_cor(label.x=.55)+
  theme_ggeffects()+
  theme(legend.position = "none")
#distplotswy<-ggarrange(irwy,dtwy,fdwy, nrow=2, ncol=2,common.legend = T )

#plots 2/8
boxwy <-ggplot(suballwy, aes(y=BRTE ,x=drought,fill=trt))+
  geom_boxplot()+
  scale_fill_viridis_d(option = "D", begin = .1, end = 1, alpha = 0.7)+
  #facet_wrap(~year,scales="fixed")+
  #geom_hline(yintercept =0,col="black")+
  labs(y="relative cover BRTE", fill="seed trt")+
  theme_ggeffects()

ggarrange(irwy,dtwy,fdwy, boxwy,nrow=2, ncol=2,common.legend = F)





mod1 <- glmmTMB(invcov.rel ~ trt * distir + (1|structure)+(1|drought),
                ziformula= ~ trt * distir + (1|structure)+(1|drought),
                data = suballca, family=beta_family())
summary(mod1)

#fesper test
mod1 <- glmmTMB(FESPER ~ drought * distir * native.cover + (1|structure),
                ziformula= ~ drought * distir *native.cover + (1|structure),
                data = suballca, family=beta_family())
summary(mod1)

### testing assumption/convergence of B-I mods
library(DHARMa)
res <- simulateResiduals(mod2)
plot(res)

##without extra 0's 
suballca2 <- suballca %>% filter(fesper.seeded=="1")
mod1.2 <- glmmTMB(invcov.rel ~ trt * drought * distir * native.cover + (1|structure),
                ziformula= ~ trt * drought * distir *native.cover + (1|structure),
                data = suballca2, family=beta_family())
mod1.2 <- glmmTMB(invcov.rel ~ trt * distir + (1|structure)+(1|drought),
                ziformula= ~ trt * distir + (1|structure)+(1|drought),
                data = suballca2, family=beta_family())
summary(mod1.2)

invloc <- read.csv("data/invasion_loc23.csv")
# suballwy2 <- merge(suballwy, invloc
# suballwy2 <- suballwy2 %>% filter(invaded=="1")
mod2 <- glmmTMB(BRTE ~ trt * drought * distir* nativecov.plot + (1|block),
                ziformula= ~ trt * drought * distir* nativecov.plot + (1|block), #does not run with distir???
                data = suballwy, family=beta_family())
summary(mod2)
