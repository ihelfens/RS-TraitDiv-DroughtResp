# FILE: R_resp_Fdiv_def_rec_pub.R
# DESCRIPTION: This script performs linear modeling on remotely-sensed trait-based functional diversity data. 
# We indicated the chosen model, tested the influence of the order of the contributions and prepared the data for further analyses and creating figures.
# Please find more information on the data acquisition and the analysis in the corresponding study. 

# USAGE: Load the script in R and ensure the CSV file is in the same directory.
#        The CSV file is named 'R_resp_Fdiv_eve_reg.csv'.

# Authors: I.S. Helfenstein, B.S. Schmid
# Date: 2024

###################################################################

# Load libraries
library(ggplot2) # for plotting
library(ggpubr) # for plotting
library(tidyverse) # for data manipulation
library(broom) # for tidying up model output
library(AICcmodavg) # for model selection
library(dplyr) # for data manipulation
library(reshape2) # for data manipulation
library(rgl) # 3D scatterplot
library(hexbin) # for hexbin plot
library(rstudioapi) 

# prepare workspace
options(digits=4) # set number of post-decimal digits in output to 4
rm(list=ls()) # clear workspace

# set working directory(replace as needed)
setwd("C:/Users/username/foldername")

## Load and preprocess data
# Richness, Evenness binned to 20 bins and 21 Regions (#8400 bins)
datfile <- read.csv("R_resp_Fdiv_reg.csv")
summary(datfile)
dat <- na.omit(datfile)

# Drought response
dat$rst <- dat$meanResistance
dat$rcv <- dat$meanRecovery
dat$rsl <- dat$meanResilience

# Diversity
dat$ric <- dat$meanRic
dat$RIC <- factor(dat$RicBin);nlevels(dat$RIC)
# Steps with the bins. ric is mean richness and the factor is the steps.
dat$eve <- dat$meanEve
dat$EVE <- factor(dat$EveBin);nlevels(dat$EVE)

# Plot steps
plot(y=dat$ric,x=dat$RicBin)
plot(y=dat$eve,x=dat$EveBin)

# Transformations
dat$logric <- log(dat$ric)
dat$ric2 <- dat$ric^2 # 
dat$ric3 <- dat$ric^3

dat$logeve <- log(dat$eve)
dat$eve2 <- dat$eve^2
dat$eve3 <- dat$eve^3

# Regions
dat$REG <- factor(dat$Region);nlevels(dat$REG)

#######################################################################
## 1 Resistance rst
# Plot
par(mfrow=c(1,3)) # three images side by side
hist(dat$rst,n=20) # extreme values (important at the top and bottom)
plot(y=dat$rst,x=dat$ric)
plot(y=dat$rst,x=dat$eve)

# A Overview
fit1 <- lm(terms(rst~logric+RIC+(eve+eve2)+EVE+logric:(eve+eve2) 
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG
           ,keep.order=T)
           ,weight=N,data=dat) 
anova(fit1) 
par(mfrow=c(2,2))
plot(fit1)
AIC(fit1)

# B Excluding extreme values of rst
upper <- mean(dat$rst,na.rm=T)+sd(dat$rst,na.rm=T)*3
lower <- mean(dat$rst,na.rm=T)-sd(dat$rst,na.rm=T)*3
dat$rstNA <- dat$rst
is.na(dat$rstNA) <- (dat$rst<lower)|(dat$rst>upper)
par(mfrow=c(1,3))
hist(dat$rstNA,n=20)
plot(y=dat$rstNA,x=dat$ric)
plot(y=dat$rstNA,x=dat$eve)

fit2 <- lm(terms(rstNA~logric+RIC+(eve+eve2)+EVE+logric:(eve+eve2)
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit2)
par(mfrow=c(2,2))
plot(fit2) 
AIC(fit2)

# C Stepwise simplified model
fit3 <- lm(terms(rst~logric+(eve+eve2)+logric:(eve+eve2)
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit3)
AIC(fit3)

fit4 <- lm(terms(rst~logric+(eve+eve2)+logric:(eve+eve2)
                 +REG+logric:REG+(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit4)
AIC(fit4)

fit5 <- lm(terms(rst~logric+(eve+eve2)+logric:(eve+eve2)
                 +REG+logric:REG+(eve):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit5)
AIC(fit5)

fit6 <- lm(terms(rst~logric+(eve+eve2)
                 +REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6)
AIC(fit6)

# D final model (i)
fitrst <- lm(terms((rst)~logric+(eve+eve2)
                 +REG+logric:REG+eve:REG+eve2:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fitrst)
AIC(fitrst)

summary(fitrst)
# This is the final model (Figure 6)

# more tests
# Different fitting sequences to see how this changes SS
fit7 <- lm(terms(rst~(eve+eve2)+logric
                 +REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit7)

fit8 <- lm(terms(rst~REG
                 +logric+(eve+eve2)
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit8)
par(mfrow=c(2,2))
plot(fit8)

# E Normalized for Regions (Figure 5)
fittemp <- lm(terms(rst~REG 
                    ,keep.order=T)
              ,weight=N,data=dat)
yadj <- residuals(fittemp)+mean(dat$rst)
#adjusted y variable, residual + overall mean
fitted8a <- lm(terms(yadj~logric+(eve+eve2)
                     ,keep.order=T)
               ,weight=N,data=dat)
anova(fitted8a)
plot3d(z=fitted(fitted8a),x=log(dat$ric),y=dat$eve)

# Use region as random effect
library(nlme) 
vf1 <- varFunc(~1/N)
fit9 <- lme(rst~logric+(eve+eve2),random=~1|REG,weights=vf1,na.action=na.omit,data=dat)

fitted9 <- fitted(fit9)
plot3d(z=dat$fitted9,x=dat$ric,y=dat$eve)

# F Fit per region (Figure S6 and S7)
fit10 <- lm(terms(rst~REG
                 +logric+REG:logric
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit10)

plot(x=dat$ric,y=fitted(fit10))

fit11 <- lm(terms(rst~REG
                  +(eve+eve2)
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit11)

plot(x=dat$eve,y=fitted(fit11))

#######################################################################
## 2 Recovery rcv
# Plot
par(mfrow=c(1,3))
hist(dat$rcv,n=20)
plot(y=dat$rcv,x=dat$ric)
plot(y=dat$rcv,x=dat$eve)

# A Overview
fit1 <- lm(terms(rcv~logric+RIC+(eve+eve2)+EVE+logric:(eve+eve2) 
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG 
                 ,keep.order=T)
           ,weight=N,data=dat) 
anova(fit1)  
par(mfrow=c(2,2))
plot(fit1) 
AIC(fit1)

# Interactions
fit5 <- lm(terms(rcv~logric+(eve+eve2)+logric:(eve+eve2)
                 +REG+logric:REG+eve:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit5)
AIC(fit5)

fit6 <- lm(terms(rcv~logric+(eve+eve2)
                 +REG+logric:REG+eve:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6)
AIC(fit6)

# D final model (ii)
fitrcv <- lm(terms(rcv~logric+(eve+eve2)
                  +REG+logric:REG+eve:REG+eve2:REG
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fitrcv)
AIC(fitrcv)

summary(fitrcv)
# This is the final model (Figure 6)

# more tests
# Different fitting sequences to see how this changes SS
fit7 <- lm(terms(rcv~(eve+eve2)+logric
                 +REG+logric:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit7)

fit8 <- lm(terms(rcv~REG
                 +logric+(eve+eve2)
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit8)
par(mfrow=c(2,2))
plot(fit8)

# E Normalized for Regions (Figure 5)
fittemp <- lm(terms(rcv~REG 
                    ,keep.order=T)
              ,weight=N,data=dat)
yadj <- residuals(fittemp)+mean(dat$rst)
fitted8a <- lm(terms(yadj~logric+(eve+eve2)
                     ,keep.order=T)
               ,weight=N,data=dat)
anova(fitted8a)
plot3d(z=fitted(fitted8a),x=log(dat$ric),y=dat$eve)

library(nlme) # for mixed-effects models
vf1 <- varFunc(~1/N) #weighting variable as variable function; gewicht ist 1/N
fit9 <- lme(rcv~logric+(eve+eve2),random=~1|REG,weights=vf1,na.action=na.omit,data=dat) 
anova(fit9)

fitted9 <- fitted(fit9)
plot3d(z=dat$fitted9,x=dat$ric,y=dat$eve)

# F Fit per region (Figure S6 and S7)
fit10 <- lm(terms(rcv~REG
                  +logric+REG:logric
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit10)

plot(x=dat$ric,y=fitted(fit10))

# Final fit Evenness
fit11 <- lm(terms(rcv~REG
                  +(eve+eve2)
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit11)

plot(x=dat$eve,y=fitted(fit11))

###########################################################################
## 3 Resilience rsl
# Plot
par(mfrow=c(1,3))
hist(dat$rsl,n=20) 
plot(y=dat$rsl,x=dat$ric)
plot(y=dat$rsl,x=dat$eve)

# A Overview
fit1 <- lm(terms(rsl~ric+logric+RIC+(eve+eve2)+EVE+ric:(eve+eve2) 
                 +REG+ric:REG+(eve+eve2):REG+ric:(eve+eve2):REG 
                 ,keep.order=T)
           ,weight=N,data=dat) 
anova(fit1) 
par(mfrow=c(2,2))
plot(fit1)
AIC(fit1)

# C Stepwise simplified model
fit3 <- lm(terms(rsl~ric+(eve+eve2)+logric:(eve+eve2)
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit3)
AIC(fit3)

fit4 <- lm(terms(rsl~ric+(eve+eve2)+ric:(eve+eve2)
                 +REG+ric:REG+(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit4)
AIC(fit4)

# Interactions
fit5 <- lm(terms(rsl~ric+(eve+eve2)+ric:(eve+eve2)
                 +REG+ric:REG+(eve):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit5)
AIC(fit5)

fit6 <- lm(terms(rsl~ric+(eve+eve2)
                 +REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6)
AIC(fit6)

# D final model (iii)
fitrsl <- lm(terms(rsl~ric+eve
                  +REG+ric:REG+eve:REG
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fitrsl)
AIC(fitrsl)

summary(fitrsl)
# This is the final model (Figure 6)

# more tests
# Different fitting sequences to see how this changes SS (i.e. variance explained)
fit7 <- lm(terms(rsl~eve+(ric+ric2)
                 +REG+eve:REG+ric:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit7)

fit8 <- lm(terms(rsl~REG
                 +(ric+ric2)+eve
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit8) 
par(mfrow=c(2,2))
plot(fit8) 

# E Normalized for Regions
fittemp <- lm(terms(rsl~REG
                    ,keep.order=T)
              ,weight=N,data=dat)
yadj <- residuals(fittemp)+mean(dat$rsl) #adjusted y variable, residual + overall mean
fitted8a <- lm(terms(yadj~ric+eve
                     ,keep.order=T)
               ,weight=N,data=dat)
anova(fitted8a)
plot3d(z=fitted(fitted8a),x=dat$ric,y=dat$eve)

# F Fit per region (Figure S6 and S7)
fit10 <- lm(terms(rsl~REG
                  +ric+REG:ric
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit10)

plot(x=dat$ric,y=fitted(fit10))
fit11 <- lm(terms(rsl~REG
                  +eve+REG:eve
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit11)

plot(x=dat$eve,y=fitted(fit11))