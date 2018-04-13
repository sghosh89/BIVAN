library(copula)
library(VineCopula)
#---------------------------------------------------------------------------------------
# This function GetNoise generates noises with a copula structure as you specified
# Input :
#        N : number of points drawn for a copula (C,G,J,F,SC,SG,SJ)
#        fcode : familycode of the desired copula [within c(3:6,13,14,16)]
#        corcoef : a number which may be Kendall's Tau or Spearman's Rho
#        method : a character of spearman or kendall
#        ploton : logical to genarate an optional plot
# Output :
#       list of two : a N by 2 noise matrix , it's qnorm transformed form and parameter of the copula
#---------------------------------------------------------------------------------------
GetNoise<-function(N,fcode,corcoef,method,ploton){

  if (fcode %in% c(3,13)){
    tgcop<-claytonCopula(3,2)
  }else if(fcode %in% c(4,14)){
    tgcop<-gumbelCopula(3,2)
  }else if(fcode %in% c(6,16)){
    tgcop<-joeCopula(3,2)
  }else if(fcode==5){
    tgcop<-frankCopula(3,2)
  }else{
    stop("fcode is not in c(3:6,13,14,16)")
  }
  
  if(method=="spearman" && fcode %in% c(3:5,13,14)){
    param<-iRho(copula = tgcop, rho = corcoef)
  }else if(method=="spearman" && !fcode %in% c(3:5,13,14)){
    warning("fcode not compatible with spearman iRho",immediate.=T,call.=T)
  }else if(method=="kendall"){
    param<-iTau(copula = tgcop, tau = corcoef)
  }else{
    param<-NA
    warning("specify method",immediate.=T,call.=T)
  }
  
  noisecop<-BiCopSim(N=N, family=fcode, par=param)
  
  if(ploton==T){
    plot(noisecop[,1],noisecop[,2],col="blue")
  }

# apply qnorm on noisecop to get normal distribution of each marginal
  noise_q<-qnorm(noisecop)
  
  return(list(noise=noisecop,noise_q=noise_q,param=param))
}

# Check the function
#s<-GetNoise(N=100,fcode=3,corcoef=0.5,method="kendall",ploton=T)

#s1<-s$noise_q
#hist(s1[,1],breaks=1000) # check if normal?
#hist(s1[,2],breaks=1000) # check if normal?

#tgcop<-claytonCopula(3,2)
#tgcop2<-rotCopula(tgcop)
#iRho(tgcop,0.5)==iRho(tgcop2,0.5) # check if it is True?
#iTau(tgcop,0.5)==iTau(tgcop2,0.5) # check if it is True?
#----------------------------------------------------------------------------------------
#Simulates a model for use in understanding Moran effects on population structure
#
#Args
#cons         An autocorrelation coefficient for the model (one number, |cons|<1)
#p0           Initial conditions - length 2 vector, default c(0,0)
#noise       An N by 2 matrix of environment variables in two habitat patches through time 
#                             (this should be noise_q from output of GetNoise function)
#burnin       The number of population times to throw away to get rid of transient behavior in the
#               dynamics. Default 500.

#Output
#An N+1 by 2 matrix of populations through time, p0 is the first row

Simulator_Cause4copula<-function(cons,p0=c(0,0),noise,burnin=500){
  
  N<-dim(noise)[1]
  cons2<-sqrt(1-cons^2)
  res<-matrix(NA,N+1,2)
  res[1,]<-p0
  
  for (counter in 2:(N+1)){
    res[counter,]<-(res[counter-1,]*cons)+(noise[counter-1,]*cons2)
  }
  
  if (burnin>0){
    res<-res[-(1:burnin),]
  }
  
  return(res)
  
} 
#-----------------------------------------------------------------------------------------------
#s<-GetNoise(N=10000,fcode=3,corcoef=0.8,method="spearman",ploton=T)
#s2<-Simulator_Cause4copula(cons=0.5,p0=c(0,0),noise=s$noise_q,burnin=1000)
#hist(s2[,1],breaks=1000) # it's normal
#s2_p<-pnorm(s2)
#hist(s2_p[,1],breaks=1000) # it's uniform
#points(s2_p[,1],s2_p[,2],col="red")

# Now compare between two correlations coming from noise_q and s2_p
#cor.test(s$noise_q[,1],s$noise_q[,2],method = "spearman")$estimate
#cor.test(s2_p[,1],s2_p[,2],method = "spearman")$estimate

#cor.test(s$noise_q[,1],s$noise_q[,2],method = "kendall")$estimate
#cor.test(s2_p[,1],s2_p[,2],method = "kendall")$estimate

#cor.test(s$noise_q[,1],s$noise_q[,2],method = "pearson")$estimate
#cor.test(s2_p[,1],s2_p[,2],method = "pearson")$estimate





