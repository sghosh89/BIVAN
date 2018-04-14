library(copula)
library(VineCopula)
source("CopulaFunctions_flexible.R")
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
  
  return(list(noise_c=noisecop,noise_q=noise_q,param=param))
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
#               dynamics. Default 1000.

#Output
#An N+1 by 2 matrix of populations through time, p0 is the first row in copula space : pop_c
# pop_q : this is similar but with normal distribution marginal

Simulator_Cause4copula<-function(cons,p0=c(0,0),noise,burnin=1000){
  
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
  
  res2<-pnorm(res) # convert into copula space
  return(list(pop_c=res2,pop_q=res))
  
} 
#-----------------------------------------------------------------------------------------------
#----------------------------------------------------------
# This function gives mean,lowCI,upCI of a vector
#
MCI<-function(x){
  m<-mean(x)
  se<-sd(x)/sqrt(length(x))   
  return(c(m-1.96*se,m,m+1.96*se))
}
#--------------------------------------------------------------------------------
#function to get a comparison table 
comp<-function(s,s2){
  
  comp<-matrix(NA,nrow=3,ncol=2)
  rownames(comp)<-c("spearman","kendall","pearson")
  colnames(comp)<-c("cor_noise","cor_pop")
  # Now compare between two correlations coming from noise_q and s2_p
  comp[1,1]<-cor(s$noise_c[,1],s$noise_c[,2],method = "spearman")
  comp[1,2]<-cor(s2$pop_c[,1],s2$pop_c[,2],method = "spearman")
  
  comp[2,1]<-cor(s$noise_c[,1],s$noise_c[,2],method = "kendall")
  comp[2,2]<-cor(s2$pop_c[,1],s2$pop_c[,2],method = "kendall")
  
  comp[3,1]<-cor(s$noise_q[,1],s$noise_q[,2],method = "pearson")
  comp[3,2]<-cor(s2$pop_q[,1],s2$pop_q[,2],method = "pearson")
  
  return(comp)
}

# when noise comes from a clayton cop with spearmancor=0.8
#s<-GetNoise(N=10000,fcode=3,corcoef=0.8,method="spearman",ploton=T)
#s2<-Simulator_Cause4copula(cons=0.5,p0=c(0,0),noise=s$noise_q,burnin=1000)
#hist(s2$pop_q[,1],breaks=1000) # it's normal
#hist(s2$pop_c[,1],breaks=1000) # it's uniform
#points(s2$pop_c[,1],s2$pop_c[,2],col="red")
#comp(s=s,s2=s2)

#
# Input :
#       numsim : a number over which desired stat (Spearman, Kendall, Pearson) called for 
#       fcode : family of copula from where noise is genarated initially [within c(3:6,13,14,16)]
#       method : a character : either "spearman" or "kendall"
#       lb : lower bound for Non-parametric stat function (Default=0)
#       ub : upper bound for Non-parametric stat function (Default =0.1)

Plotter_Cause4copula_stat<-function(numsim,fcode,method,lb=0,ub=0.1){
  
  corcoef_list<-seq(from=0.1,to=0.9,by=0.1)
  
  # initialize
  S_noise_mat<-matrix(NA,nrow=length(corcoef_list),ncol=3) # Spearman Correlation matrix for noise
  colnames(S_noise_mat)<-c("lowCI","mean","upCI")  
  K_noise_mat<-S_noise_mat  # Kendall Correlation matrix
  P_noise_mat<-S_noise_mat  # Pearson correlation matrix
  S_pop_mat<-S_noise_mat    # Spearman Correlation matrix for population
  K_pop_mat<-S_noise_mat
  P_pop_mat<-S_noise_mat
  Corl_noise_mat<-S_noise_mat # Corl stat matrix for noise
  Corl_pop_mat<-S_noise_mat
  Coru_noise_mat<-S_noise_mat # Coru stat matrix for noise
  Coru_pop_mat<-S_noise_mat
  Pl_noise_mat<-S_noise_mat 
  Pl_pop_mat<-S_noise_mat    # Pl stat matrix for population
  Pu_noise_mat<-S_noise_mat
  Pu_pop_mat<-S_noise_mat
  D2u_noise_mat<-S_noise_mat # D2u stat matrix for noise
  D2u_pop_mat<-S_noise_mat
  D2l_noise_mat<-S_noise_mat # D2u stat matrix for noise
  D2l_pop_mat<-S_noise_mat
  
  pval_S<-c() # an empty vector to store p values from t test of Spearman Cor. from noise and population
  pval_K<-c()
  pval_P<-c()
  pval_Corl<-c()
  pval_Coru<-c()
  pval_Pl<-c()
  pval_Pu<-c()
  pval_D2u<-c()
  pval_D2l<-c()
  
  for(counter in c(1:length(corcoef_list))){
    corcoef<-corcoef_list[counter]
    #cat("corcoef=",corcoef,"\n")
    S_noise<-c()
    K_noise<-c()
    P_noise<-c()
    S_pop<-c()
    K_pop<-c()
    P_pop<-c()
    Corl_noise<-c()
    Corl_pop<-c()
    Coru_noise<-c()
    Coru_pop<-c()
    Pl_noise<-c()
    Pl_pop<-c()
    Pu_noise<-c()
    Pu_pop<-c()
    D2u_noise<-c()
    D2u_pop<-c()
    D2l_noise<-c()
    D2l_pop<-c()
    for(i in 1:numsim){
      #cat("i=",i,"\n")
      s<-GetNoise(N=10000,fcode=fcode,corcoef=corcoef,method=method,ploton=F)
      s2<-Simulator_Cause4copula(cons=0.5,p0=c(0,0),noise=s$noise_q,burnin=1000)
      z<-as.data.frame(comp(s=s,s2=s2))
      S_noise<-c(S_noise,z$cor_noise[1])
      K_noise<-c(K_noise,z$cor_noise[2])
      P_noise<-c(P_noise,z$cor_noise[3])
      S_pop<-c(S_pop,z$cor_pop[1])
      K_pop<-c(K_pop,z$cor_pop[2])
      P_pop<-c(P_pop,z$cor_pop[3])
      
      Corl_noise<-c(Corl_noise,Corbds(vi = s$noise_c[,1],vj = s$noise_c[,2],lb = lb,ub = ub))
      Corl_pop<-c(Corl_pop,Corbds(vi = s2$pop_c[,1],vj = s2$pop_c[,2],lb = lb,ub = ub))
      
      Coru_noise<-c(Coru_noise,Corbds(vi = s$noise_c[,1],vj = s$noise_c[,2],lb = 1-ub,ub = 1-lb))
      Coru_pop<-c(Coru_pop,Corbds(vi = s2$pop_c[,1],vj = s2$pop_c[,2],lb = 1-ub,ub = 1-lb))
      
      Pl_noise<-c(Pl_noise,Pbds(vi = s$noise_c[,1],vj = s$noise_c[,2],lb = lb,ub = ub)$abs_res)
      Pl_pop<-c(Pl_pop,Pbds(vi = s2$pop_c[,1],vj = s2$pop_c[,2],lb = lb,ub = ub)$abs_res)
      
      Pu_noise<-c(Pu_noise,Pbds(vi = s$noise_c[,1],vj = s$noise_c[,2],lb = 1-ub,ub = 1-lb)$abs_res)
      Pu_pop<-c(Pu_pop,Pbds(vi = s2$pop_c[,1],vj = s2$pop_c[,2],lb = 1-ub,ub = 1-lb)$abs_res)
      
      D2u_noise<-c(D2u_noise,D2bds(vi = s$noise_c[,1],vj = s$noise_c[,2],lb = 1-ub,ub = 1-lb))
      D2u_pop<-c(D2u_pop,D2bds(vi = s2$pop_c[,1],vj = s2$pop_c[,2],lb = 1-ub,ub = 1-lb))
      
      D2l_noise<-c(D2l_noise,D2bds(vi = s$noise_c[,1],vj = s$noise_c[,2],lb = lb,ub = ub))
      D2l_pop<-c(D2l_pop,D2bds(vi = s2$pop_c[,1],vj = s2$pop_c[,2],lb = lb,ub = ub))
      
    }
    S_noise_mat[counter,]<-MCI(S_noise)
    K_noise_mat[counter,]<-MCI(K_noise)
    P_noise_mat[counter,]<-MCI(P_noise)
    S_pop_mat[counter,]<-MCI(S_pop)
    K_pop_mat[counter,]<-MCI(K_pop)
    P_pop_mat[counter,]<-MCI(P_pop)
    Corl_noise_mat[counter,]<-MCI(Corl_noise)
    Corl_pop_mat[counter,]<-MCI(Corl_pop)
    Coru_noise_mat[counter,]<-MCI(Coru_noise)
    Coru_pop_mat[counter,]<-MCI(Coru_pop)
    Pl_noise_mat[counter,]<-MCI(Pl_noise)
    Pl_pop_mat[counter,]<-MCI(Pl_pop)
    Pu_noise_mat[counter,]<-MCI(Pu_noise)
    Pu_pop_mat[counter,]<-MCI(Pu_pop)
    D2u_noise_mat[counter,]<-MCI(D2u_noise)
    D2u_pop_mat[counter,]<-MCI(D2u_pop)
    D2l_noise_mat[counter,]<-MCI(D2l_noise)
    D2l_pop_mat[counter,]<-MCI(D2l_pop)
    
    pval_S<-c(pval_S,t.test(S_noise,S_pop,alternative="two.sided",paired=T)$p.value)
    pval_K<-c(pval_K,t.test(K_noise,K_pop,alternative="two.sided",paired=T)$p.value)
    pval_P<-c(pval_P,t.test(P_noise,P_pop,alternative="two.sided",paired=T)$p.value)
    pval_Corl<-c(pval_Corl,t.test(Corl_noise,Corl_pop,alternative="two.sided",paired=T)$p.value)
    pval_Coru<-c(pval_Coru,t.test(Coru_noise,Coru_pop,alternative="two.sided",paired=T)$p.value)
    pval_Pl<-c(pval_Pl,t.test(Pl_noise,Pl_pop,alternative="two.sided",paired=T)$p.value)
    pval_Pu<-c(pval_Pu,t.test(Pu_noise,Pu_pop,alternative="two.sided",paired=T)$p.value)
    pval_D2u<-c(pval_D2u,t.test(D2u_noise,D2u_pop,alternative="two.sided",paired=T)$p.value)
    pval_D2l<-c(pval_D2l,t.test(D2l_noise,D2l_pop,alternative="two.sided",paired=T)$p.value)
    
  }
  
  if(method=="spearman"){
    xlabel<-"Spearman's Rho"
  }else if(method=="kendall"){
    xlabel<-"Kendall's Tau"
  }else{
    warning("specify method",immediate.=T,call.=T)
  }
  
  op<-par(mfrow=c(3,3),mar=c(3,3.5,3,3.5), mgp=c(1.5,0.5,0))
  plot(corcoef_list,S_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Spearman",xlim=c(0,1),ylim=c(0,1))
  arrows(corcoef_list,S_noise_mat[,1],corcoef_list,S_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,S_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,S_pop_mat[,1],corcoef_list,S_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  points(corcoef_list,pval_S,col="purple")
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  #  legend("topleft", c("noise","population"), fill=c('red', 'blue'), horiz=T, bty='n')
  
  plot(corcoef_list,K_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Kendall",xlim=c(0,1),ylim=c(0,1))
  arrows(corcoef_list,K_noise_mat[,1],corcoef_list,K_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,K_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,K_pop_mat[,1],corcoef_list,K_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  points(corcoef_list,pval_K,col="purple")
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,P_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Pearson",xlim=c(0,1),ylim=c(0,1))
  arrows(corcoef_list,P_noise_mat[,1],corcoef_list,P_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,P_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,P_pop_mat[,1],corcoef_list,P_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  points(corcoef_list,pval_P,col="purple")
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,Corl_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Corl",
       xlim=c(0,1),ylim=c(0,0.1+max(Corl_noise_mat[,2],Corl_pop_mat[,2])))
  arrows(corcoef_list,Corl_noise_mat[,1],corcoef_list,Corl_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Corl_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,Corl_pop_mat[,1],corcoef_list,Corl_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Corl,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,Coru_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Coru",
       xlim=c(0,1),ylim=c(0,0.1+max(Coru_noise_mat[,2],Coru_pop_mat[,2])))
  arrows(corcoef_list,Coru_noise_mat[,1],corcoef_list,Coru_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Coru_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,Coru_pop_mat[,1],corcoef_list,Coru_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Coru,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4, col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,Pl_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Pl",
       xlim=c(0,1),ylim=c(0,0.02+max(Pl_noise_mat[,2],Pl_pop_mat[,2])))
  arrows(corcoef_list,Pl_noise_mat[,1],corcoef_list,Pl_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Pl_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,Pl_pop_mat[,1],corcoef_list,Pl_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Pl,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,Pu_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Pu",
       xlim=c(0,1),ylim=c(0,0.02+max(Pu_noise_mat[,2],Pu_pop_mat[,2])))
  arrows(corcoef_list,Pu_noise_mat[,1],corcoef_list,Pu_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Pu_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,Pu_pop_mat[,1],corcoef_list,Pu_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Pu,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,D2u_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="D2u",
       xlim=c(0,1),ylim=c(0,0.002+max(D2u_noise_mat[,2],D2u_pop_mat[,2])))
  arrows(corcoef_list,D2u_noise_mat[,1],corcoef_list,D2u_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,D2u_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,D2u_pop_mat[,1],corcoef_list,D2u_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_D2u,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  plot(corcoef_list,D2l_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="D2l",
       xlim=c(0,1),ylim=c(0,0.002+max(D2l_noise_mat[,2],D2l_pop_mat[,2])))
  arrows(corcoef_list,D2l_noise_mat[,1],corcoef_list,D2l_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,D2l_pop_mat[,2],cex=0.5,col="blue")
  arrows(corcoef_list,D2l_pop_mat[,1],corcoef_list,D2l_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_D2l,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  #pval_range<-c(0,1)
  #print(pretty(pval_range))
  #axis(side=4, at = pretty(0:1),col='purple',col.axis="purple")
  axis(side=4,col='purple',col.axis="purple")
  mtext(side = 4, line = 1.5, 'p values', col='purple')
  
  par(op)
  
  op2<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("top", c("noise","population"), col = c("red", "blue"),
         cex = 0.8,lwd = 1, lty = 1, xpd = TRUE, horiz = TRUE, inset = c(0,0), 
         bty = "n") 
  par(op2)
  
}

