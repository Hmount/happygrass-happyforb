rain <- read.csv("data/site_data/rainGaugeTrt.csv", header=T)
rain$rain_mm_1 <- rain$rain_in_2021_5_31*25.4
rain$rain_mm_2 <- rain$rain_in_2021_6_9*25.4

library(lme4)
library(lmerTest)

m1 <- lmer(rain_mm_1 ~ trt + (1|block), data=rain)
summary(m1)
boxplot(rain_mm_1 ~ trt, data=rain, ylab="mm rain", xlab="Treatment")
t.test(rain_mm_1 ~ trt, data=rain)

#calculate percent reduction = 64%
7.636933 - 2.726267
(4.910666/7.636933)*100