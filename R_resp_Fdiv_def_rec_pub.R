# FILE: R_resp_Fdiv_def_rec_pub.R
# DESCRIPTION: This script performs various linear and mixed-effects model analyses. 
# Please find more information on the data acquisition and the analysis in the corresponding study. 
# It also exports the results to text files for further analysis and figure creation in MATLAB.

# USAGE: Load the script in R and ensure the CSV file is in the same directory.
#        The CSV file is named 'R_resp_Fdiv_eve_reg.csv'.

# Author: I.S. Helfenstein
# Date: 2024

###################################################################
# Load necessary libraries
library(car)
library(rgl)
library(nlme)

options(digits=4) # set number of post-decimal digits in output to 4
rm(list=ls()) # clear workspace

# set working directory
setwd("C:/Users/ihelfens/Documents/RS_Resilience/data/export_data")

# Load and preprocess data
datfile <- read.csv("R_resp_Fdiv_eve_reg.csv")
dat <- na.omit(datfile)
summary(dat) # Richness, Evenness data binned to 20 bins and 21 Regions (#8400)

# Define variables
dat$rst <- dat$meanResistance
dat$rcv <- dat$meanRecovery
dat$rsl <- dat$meanResilience
dat$ric <- dat$meanRic
dat$RIC <- factor(dat$RicBin); nlevels(dat$RIC) # Number of levels in RicBin
dat$eve <- dat$meanEve
dat$EVE <- factor(dat$EveBin); nlevels(dat$EVE) # Number of levels in EveBin
dat$REG <- factor(dat$Region); nlevels(dat$REG) # Number of levels in Region

# Plot richness and evenness against their bins
plot(y=dat$ric, x=dat$RicBin)
plot(y=dat$eve, x=dat$EveBin)

# Transform variables
dat$logric <- log(dat$ric) # Log-transformed richness
dat$ric2 <- dat$ric^2 # Squared richness
dat$ric3 <- dat$ric^3 # Cubed richness
dat$logeve <- log(dat$eve) # Log-transformed evenness
dat$eve2 <- dat$eve^2 # Squared evenness
dat$eve3 <- dat$eve^3 # Cubed evenness

##########################################################
# RESISTANCE rst
# Drought response (Resistance) ~ richness + evenness + richness:evenness
par(mfrow=c(1,3)) # Arrange three plots in one row
hist(dat$rst, n=20) # Histogram of resistance with 20 bins
plot(y=dat$rst, x=dat$ric) # Plot resistance against richness
plot(y=dat$rst, x=dat$eve) # Plot resistance against evenness

fit1 <- lm(terms(rst~logric+RIC+(eve+eve2)+EVE+logric:(eve+eve2)
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit1)
par(mfrow=c(2,2))
plot(fit1)
AIC(fit1)

# Excluding extreme values of rst: very little effect
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

# Simplify model: Remove Ric and Eve factors (bins)
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

# Completely without interactions
fit6 <- lm(terms(rst~logric+(eve+eve2)
                 +REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6)
AIC(fit6)

# Chosen model
fit6A <- lm(terms((rst)~logric+(eve+eve2)
                 +REG+logric:REG+eve:REG+eve2:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6A)
par(mfrow=c(2,2))
plot(fit6A)
AIC(fit6A)

summary(fit6A) # Summary of the model, r2 = 0.865

# Different fitting sequences
fit7 <- lm(terms(rst~(eve+eve2)+logric
                 +REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit7)

# Different fitting sequences
fit8 <- lm(terms(rst~REG
                 +logric+(eve+eve2)
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit8)
par(mfrow=c(2,2))
plot(fit8)

# Only regions for illustration of fit
fittemp <- lm(terms(rst~REG
                    ,keep.order=T)
              ,weight=N,data=dat)
yadj <- residuals(fittemp)+mean(dat$rst)
fitted8a <- lm(terms(yadj~logric+(eve+eve2)
                     ,keep.order=T)
               ,weight=N,data=dat)
anova(fitted8a)

# Illustrations
plot3d(z=fitted(fitted8a),x=log(dat$ric),y=dat$eve)

# Use region as random effect
vf1 <- varFunc(~1/N)
fit9 <- lme(rst~logric+(eve+eve2),random=~1|REG,weights=vf1,na.action=na.omit,data=dat)
anova(fit9)

fitted9 <- fitted(fit9)
plot3d(z=fitted9,x=dat$ric,y=dat$eve)

# Plot fit per region Richness
fit10 <- lm(terms(rst~REG
                 +logric+REG:logric
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit10)

plot(x=dat$ric,y=fitted(fit10))

# Plot fit per region Evenness
fit11 <- lm(terms(rst~REG
                  +(eve+eve2)
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit11)

plot(x=dat$eve,y=fitted(fit11))

#######################################################################################################################################
# Resilience rsl
# drought response (Resistance) ~ richness + evenness + richness:evenness.
par(mfrow=c(1,3)) # drei bilder nebeneinander
hist(dat$rsl,n=20) # viele zimlich extreme werte (wichtig oben und unten!)
plot(y=dat$rsl,x=dat$ric)
plot(y=dat$rsl,x=dat$eve)

# Resilience rsl
fit1 <- lm(terms(rsl~ric+logric+RIC+(eve+eve2)+EVE+ric:(eve+eve2)
                 +REG+ric:REG+(eve+eve2):REG+ric:(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit1)
par(mfrow=c(2,2))
plot(fit1)
AIC(fit1)

# Simplify model: Remove Ric and Eve factors (bins)
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

# Without interactions
fit5 <- lm(terms(rsl~ric+(eve+eve2)+ric:(eve+eve2)
                 +REG+ric:REG+(eve):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit5)
AIC(fit5)

# Completely without interactions
fit6 <- lm(terms(rsl~ric+(eve+eve2)
                 +REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6)
AIC(fit6)

# Chosen model
fit6A <- lm(terms(rsl~ric+eve
                  +REG+ric:REG+eve:REG
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit6A)
AIC(fit6A)

# Total model explains a lot:
summary(fit6A) # R2 = 0.74

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

# Only regions for illustration of fit
fittemp <- lm(terms(rsl~REG
                    ,keep.order=T)
              ,weight=N,data=dat)
yadj <- residuals(fittemp)+mean(dat$rsl)
fitted8a <- lm(terms(yadj~ric+eve
                     ,keep.order=T)
               ,weight=N,data=dat)
anova(fitted8a)
plot3d(z=fitted(fitted8a),x=dat$ric,y=dat$eve)

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

#######################################################################################################################################
# Recovery rcv
# drought response (Recovery) ~ richness + evenness + richness:evenness.
par(mfrow=c(1,3))
hist(dat$rcv,n=20)
plot(y=dat$rcv,x=dat$ric)
plot(y=dat$rcv,x=dat$eve)

# Recovery rcv
fit1 <- lm(terms(rcv~logric+RIC+(eve+eve2)+EVE+logric:(eve+eve2)
                 +REG+logric:REG+(eve+eve2):REG+logric:(eve+eve2):REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit1)
par(mfrow=c(2,2))
plot(fit1)
AIC(fit1)

# Without interactions
fit5 <- lm(terms(rcv~logric+(eve+eve2)+logric:(eve+eve2)
                 +REG+logric:REG+eve:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit5)
AIC(fit5)

# Completely without interactions
fit6 <- lm(terms(rcv~logric+(eve+eve2)
                 +REG+logric:REG+eve:REG
                 ,keep.order=T)
           ,weight=N,data=dat)
anova(fit6)
AIC(fit6)

# Chosen model
fit6A <- lm(terms(rcv~logric+(eve+eve2)
                  +REG+logric:REG+eve:REG+eve2:REG
                  ,keep.order=T)
            ,weight=N,data=dat)
anova(fit6A)
AIC(fit6A)

# We can see that the total model explains a lot:
summary(fit6A) # R2 = 0.732

# Different fitting sequences to see how this changes SS (i.e. variance explained)
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

# Only regions for illustration of fit
fittemp <- lm(terms(rcv~REG # einfach nur region
                    ,keep.order=T)
              ,weight=N,data=dat)
yadj <- residuals(fittemp)+mean(dat$rst)
# Fit a linear model with log-transformed richness and evenness variables
fitted8a <- lm(terms(yadj~logric+(eve+eve2)
                     ,keep.order=T)
               ,weight=N,data=dat)
anova(fitted8a) # Perform ANOVA on the fitted model
plot3d(z=fitted(fitted8a),x=log(dat$ric),y=dat$eve) # Plot the fitted values in 3D

# Use region as a random effect with all environmental factors included
vf1 <- varFunc(~1/N) # Define a weighting variable as a variable function
fit9 <- lme(rcv~logric+(eve+eve2),random=~1|REG,weights=vf1,na.action=na.omit,data=dat) # Fit a mixed-effects model with region as a random effect
anova(fit9) # Perform ANOVA on the fitted mixed-effects model

fitted9 <- fitted(fit9) # Get the fitted values from the mixed-effects model
plot3d(z=dat$fitted9,x=dat$ric,y=dat$eve) # Plot the fitted values in 3D

fit10 <- lm(terms(rcv~REG
                  +logric+REG:logric
                  ,keep.order=T)
            ,weight=N,data=dat) # Fit a linear model with region and log-transformed richness variables
anova(fit10) # Perform ANOVA on the fitted model
#Region only takes little away from richness and evenness

plot(x=dat$ric,y=fitted(fit10)) # Plot the fitted values against richness

fit11 <- lm(terms(rcv~REG
                  +(eve+eve2)
                  ,keep.order=T)
            ,weight=N,data=dat) # Fit a linear model with region and evenness variables
anova(fit11) # Perform ANOVA on the fitted model
# Region only takes little away from richness and evenness

plot(x=dat$eve,y=fitted(fit11)) # Plot the fitted values against evenness
