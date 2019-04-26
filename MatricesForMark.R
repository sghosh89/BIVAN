library(mvtnorm)
library(copula)

rm(list=ls())
set.seed(101)

#extreme left tail dependent
source("ExtremeTailDep.R")
d<-retd(1000000,2,-1)
plot(d[1:1000,1],d[1:1000,2],type="p")
h<-cor(d,method="spearman") #should about equal rhoval
h
rhoval<-h[2,1]
hist(d[,1])
hist(d[,2]) #these should look normal
var(d[,1])
var(d[,2])
cor(d)
write.csv(d,file="ExtremeLeftTailDep.csv")

#extreme right tail dependent
d<-(-d)
plot(d[1:1000,1],d[1:1000,2],type="p")
h<-cor(d,method="spearman") #should equal rhoval
h
rhoval
hist(d[,1])
hist(d[,2]) #these should look normal
var(d[,1])
var(d[,2])
cor(d)
write.csv(d,file="ExtremeRightTailDep.csv")

#covariance to use
#pearval<-0.7 #covariance, Pearson
#ncop<-normalCopula(pearval,2)
#rhoval<-rho(ncop) #spearman

#bivariate normal
ncop<-normalCopula(.5,2)
nparam<-iRho(ncop,rhoval)
ncop<-normalCopula(nparam,2)
rho(ncop) #should equal rhoval
rhoval
dobs<-rCopula(1000000,ncop)
plot(dobs[1:1000,1],dobs[1:1000,2],type='p') #should look like a normal copula
d<-qnorm(dobs)
dim(d)
cor(d,method="spearman") #should about equal rhoval
rhoval
cor(d)
hist(d[,1])
hist(d[,2]) #these should look normal
write.csv(d,file="BivariateNormal.csv")

#somewhat left-tail dependent
ccop<-claytonCopula(1,2)
cparam<-iRho(ccop,rhoval)
ccop<-claytonCopula(cparam,2)
rho(ccop) #should equal rhoval
rhoval
dobs<-rCopula(1000000,ccop)
plot(dobs[1:1000,1],dobs[1:1000,2],type='p') #should look like a clayton copula
d<-qnorm(dobs)
dim(d)
cor(d,method="spearman") #should about equal rhoval
rhoval
cor(dobs,method="spearman") #should be the same
cor(d)
hist(d[,1])
hist(d[,2]) #these should look normal
write.csv(d,file="SomewhatLeftTailDep.csv")

#somewhat right-tail dependent
dobs<-(1-dobs)
plot(dobs[1:1000,1],dobs[1:1000,2],type='p') #should look like a survival clayton copula
d<-qnorm(dobs)
dim(d)
cor(d,method="spearman") #should about equal rhoval
rhoval
cor(dobs,method="spearman") #should be the same
hist(d[,1])
hist(d[,2]) #these should look normal
cor(d)
write.csv(d,file="SomewhatRightTailDep.csv")


