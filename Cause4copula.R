library(copula)
library(VineCopula)
library(e1071)
source("CopulaFunctions_flexible.R")
source("MyBiCopGofTest.R")
#---------------------------------------------------------------------------------------
# This function GetNoise generates noises with a copula structure as you specified
# Args :
#        N : number of points drawn for a copula (C,G,J,F,SC,SG,SJ)
#        fcode : familycode of the desired copula [within c(3:6,13,14,16)]
#        corcoef : a number which may be Kendall's Tau or Spearman's Rho
#        nsd : standard deviation of noise generated
#        method : a character of spearman or kendall
#        ploton : logical to genarate an optional plot
# Output :
#       list of 3 : 
#                  noise_c : a N by 2 noise matrix ,
#                  noise_q : it's qnorm transformed form 
#                  param : parameter of the copula from where that noise is generated
#---------------------------------------------------------------------------------------
GetNoise<-function(N,fcode,corcoef,nsd,method,ploton){

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
    warning("specify method compatible with copula family",immediate.=T,call.=T)
  }
  
  noisecop<-BiCopSim(N=N, family=fcode, par=param)
  
  if(ploton==T){
    plot(noisecop[,1],noisecop[,2],col="blue")
  }

# apply qnorm on noisecop to get normal distribution of each marginal
  noise_q<-qnorm(noisecop,mean=0,sd=nsd)
  
  return(list(noise_c=noisecop,noise_q=noise_q,param=param))
}
#------------------------------------------------------------------------------
# Check the function
#s<-GetNoise(N=10000,fcode=3,corcoef=0.5,nsd=0.1,method="kendall",ploton=T)

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
#params       a vector  with 
#               if model == "ar1" : an autocorrelation coefficient for the model (one number, |cons|<1)
#               if model == "ricker" : with two number, 1. r, 2. K
                            
#p0           Initial conditions - length 2 vector, default c(0,0)
#noise       An N by 2 matrix of environment variables in two habitat patches through time 
#                             (this should be noise_q from output of GetNoise function)

#Output
# list of 2 : 
#           pop_c : An N+1 by 2 matrix of populations through time, p0 is the first row in copula space
#           pop_q : this is similar but with normal distribution marginal

Simulator_Cause4copula<-function(params,p0=c(0,0),noise,model){
  
  N<-dim(noise)[1]
  res<-matrix(NA,N+1,2)
  res[1,]<-p0
  
  
    if(model=="ar1"){
      beta<-params
      alpha<-sqrt(1-beta^2)
      for (counter in 2:(N+1)){
        res[counter,]<-(beta*res[counter-1,])+(alpha*noise[counter-1,])
      }
    }else if(model=="ricker"){
      r<-params[1]
      K<-params[2]
      for (counter in 2:(N+1)){
        res[counter,]<-res[counter-1,]*exp((r*(1-(res[counter-1,]/K)))+noise[counter-1,])
      }
    }else{
      warning("model not specified",immediate.=T,call.=T)
    }
  
  # convert into copula space
  #res2<-pnorm(res) 
  res2<-VineCopula::pobs(res)
  return(list(pop_c=res2,pop_q=res))
  
} 
#-----------------------------------------------------------------------------------------------
# The following lines are examples to draw error bar without calling arrows
#plot(x=c(1:5),y=c(1:5),ylim=c(0,6))
#se<-c(1,0.5,NA,0,0.8)
#segments(x,y-se,x,y+se)
#epsilon <- 0.02
#segments(x-epsilon,y-se,x+epsilon,y-se)
#segments(x-epsilon,y+se,x+epsilon,y+se)
#----------------------------------------------------------
# This function gives mean,lowCI,upCI of a vector
#
MCI<-function(x){
  x<-x[is.finite(x)] # only considering finite values
  m<-mean(x)
  se<-sd(x)/sqrt(length(x))   
  return(c(m-1.96*se,m,m+1.96*se))
}
#------------------------------------------------------------------------------------------------------------
#function to get a list with a comparison table, chopped data for input noise and output copula
# Args:
#      s : gives you [s$noise_c, s$noise_q : each is a N by 2 noise matrix] , param
#      s2 : gives you [s2$pop_c, s2$pop_q : each is a N+1 by 2 noise matrix] 
#     num_keep_last : an integer : number of rows you want to keep from bottom 
#                     for each of s$noise_c, s$noise_q and s2$pop_c, s2$pop_q matrix 
#
comp<-function(s,s2,num_keep_last){
  
  mod_s_noise_c<-tail(s$noise_c,num_keep_last)
  mod_s_noise_q<-tail(s$noise_q,num_keep_last)
  mod_s2_pop_c<-tail(s2$pop_c,num_keep_last)
  mod_s2_pop_q<-tail(s2$pop_q,num_keep_last)
  
  comp<-matrix(NA,nrow=3,ncol=2)
  rownames(comp)<-c("spearman","kendall","pearson")
  colnames(comp)<-c("cor_noise","cor_pop")
  
  # Now compare between spearman corln btw noise_c and pop_c
  comp[1,1]<-cor(mod_s_noise_c[,1],mod_s_noise_c[,2],method = "spearman")
  comp[1,2]<-cor(mod_s2_pop_c[,1],mod_s2_pop_c[,2],method = "spearman")
  
  # Now compare between kendall corln btw noise_c and pop_c
  comp[2,1]<-cor(mod_s_noise_c[,1],mod_s_noise_c[,2],method = "kendall")
  comp[2,2]<-cor(mod_s2_pop_c[,1],mod_s2_pop_c[,2],method = "kendall")
  
  # Now compare between pearson corln btw noise_q and pop_q
  comp[3,1]<-cor(mod_s_noise_q[,1],mod_s_noise_q[,2],method = "pearson")
  comp[3,2]<-cor(mod_s2_pop_q[,1],mod_s2_pop_q[,2],method = "pearson")
  
  rownames(mod_s_noise_c)<-c()
  rownames(mod_s_noise_q)<-c()
  rownames(mod_s2_pop_c)<-c()
  rownames(mod_s2_pop_q)<-c()
  
  last_num_keep_noise<-list(last_num_keep_noise_c=mod_s_noise_c,
                            last_num_keep_noise_q=mod_s_noise_q)
  
  last_num_keep_pop<-list(last_num_keep_pop_c=mod_s2_pop_c,
                            last_num_keep_pop_q=mod_s2_pop_q)
  
  return(list(comp=as.data.frame(comp),
              last_num_keep_noise=last_num_keep_noise,
              last_num_keep_pop=last_num_keep_pop))
}
#-------------------------------------------------------------------
# when noise comes from a clayton cop with spearmancor=0.8
#s<-GetNoise(N=1000,fcode=3,corcoef=0.3,method="spearman",ploton=T)
#s2<-Simulator_Cause4copula(params=c(0.5),p0=c(0,0),noise=s$noise_q,model="ar1")
#hist(s2$pop_q[,1],breaks=1000) # it's normal
#hist(s2$pop_c[,1],breaks=1000) # it's uniform
#plot(s2$pop_c[,1],s2$pop_c[,2],col="red")

#set.seed(seed=101)
#s<-GetNoise(N=5000,fcode=3,corcoef=0.1,method="spearman",ploton=T)
#hist(s$noise_q[,1],breaks=100) #should be normal
#params<-c(0.8,100)
#p0<-c(100,100)
#noise<-s$noise_q
#noise<-matrix(0,nrow=1000,ncol=2)
#noise<-noise*0.75
#s2<-Simulator_Cause4copula(params=params,p0=p0,noise=noise,model="ricker")
#t<-c(1:nrow(s2$pop_q))
#plot(t,s2$pop_q[,1],type="l")
#hist(s2$pop_q[,1],breaks=1000) # it should be normal, but it's not???
#hist(s2$pop_c[,1],breaks=1000) # it's uniform
#plot(s2$pop_c[,1],s2$pop_c[,2],col="red")
#zz<-comp(s=s,s2=s2)
#---------------------------------------------------------

#---------------------------------------------------------------------------------------------------
# This function computes
#  Spearman correlation btw noise_c & pop_c against spearman's rho/ kendall's tau of noise
#  Kendall correlation btw noise_c & pop_c against spearman's rho/ kendall's tau of noise
#  Pearson correlation btw noise_q & pop_q against spearman's rho/ kendall's tau of noise
#  non-parametric statistics @ extreme ends (Corl,Coru,Pl,Pu,D2u,D2l) 
#                                 against spearman's rho/ kendall's tau of noise
# on both output copula and input noise copula 
# In addition, it computes the standard error associated with total number of 
# simulations for each value computed and estimates P -values 
# from a paired t-test of the null hypothesis that the distributions of noise and output 
# values have the same mean
#
# Args :
#       N : number of points drawn for a noise copula initially
#       numsim : a number over which desired stat (Spearman, Kendall, Pearson) called for (default:50)
#       fcode : family of copula from where noise is genarated initially [within c(3:6,13,14,16)]
#       nsd : standard deviation of the noise generated
#       method : a character : either "spearman" or "kendall"
#       lb : lower bound for Non-parametric stat function (Default=0)
#       ub : upper bound for Non-parametric stat function (Default =0.1)
#       num_keep_last : number of rows you want to keep from bottom 
#                     for each of s$noise_c, s$noise_q and s2$pop_c, s2$pop_q matrix
#       resloc : folder location where the plots should be saved
#       params,p0,model : these are inputs for Simulator_Cause4copula function
#-----------------------------------------------------------------------------------------------
Sim_Cause4copula_stat<-function(N,numsim=50,fcode,nsd,method,lb=0,ub=0.1,num_keep_last,resloc,params,p0,model){
  
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
  CorlmCoru_noise_mat<-S_noise_mat # Corl-Coru stat matrix for noise
  CorlmCoru_pop_mat<-S_noise_mat   # Corl-Coru stat matrix for pop
  Pl_noise_mat<-S_noise_mat 
  Pl_pop_mat<-S_noise_mat    # Pl stat matrix for population
  Pu_noise_mat<-S_noise_mat
  Pu_pop_mat<-S_noise_mat
  PlmPu_noise_mat<-S_noise_mat # Pl-Pu stat matrix for noise
  PlmPu_pop_mat<-S_noise_mat   # Pl-Pu stat matrix for pop
  D2u_noise_mat<-S_noise_mat # D2u stat matrix for noise
  D2u_pop_mat<-S_noise_mat
  D2l_noise_mat<-S_noise_mat # D2u stat matrix for noise
  D2l_pop_mat<-S_noise_mat
  D2umD2l_noise_mat<-S_noise_mat # D2u-D2l stat matrix for noise
  D2umD2l_pop_mat<-S_noise_mat
  
  pval_S<-c() # an empty vector to store p values from t test of Spearman Cor. from noise and population
  pval_K<-c()
  pval_P<-c()
  pval_Corl<-c()
  pval_Coru<-c()
  pval_Pl<-c()
  pval_Pu<-c()
  pval_D2u<-c()
  pval_D2l<-c()
  pval_CorlmCoru<-c()
  pval_PlmPu<-c()
  pval_D2umD2l<-c()
  
  for(counter in c(1:length(corcoef_list))){
    corcoef<-corcoef_list[counter]
    #cat("--------corcoef=",corcoef,"---------\n")
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
    CorlmCoru_noise<-c()
    CorlmCoru_pop<-c()
    Pl_noise<-c()
    Pl_pop<-c()
    Pu_noise<-c()
    Pu_pop<-c()
    PlmPu_noise<-c()
    PlmPu_pop<-c()
    D2u_noise<-c()
    D2u_pop<-c()
    D2l_noise<-c()
    D2l_pop<-c()
    D2umD2l_noise<-c()
    D2umD2l_pop<-c()
    
    all_numsim_pop_patch1<-c()# to store all pop_q from each numsim for 1st patch : this should be num_keep_last times numsim points
    all_numsim_pop_patch2<-c()# to store all pop_q from each numsim for 2nd patch : this should be num_keep_last times numsim points
    
    for(i in 1:numsim){
      #cat("i=",i,"\n")
      s<-GetNoise(N=N,fcode=fcode,corcoef=corcoef,nsd=nsd,method=method,ploton=F)
      s2<-Simulator_Cause4copula(params=params,p0=p0,noise=s$noise_q,model=model)
      z<-comp(s=s,s2=s2,num_keep_last = num_keep_last)
      
      all_numsim_pop_patch1<-c(all_numsim_pop_patch1,z$last_num_keep_pop$last_num_keep_pop_q[,1])
      all_numsim_pop_patch2<-c(all_numsim_pop_patch2,z$last_num_keep_pop$last_num_keep_pop_q[,2])
      
      S_noise<-c(S_noise,z$comp$cor_noise[1])
      K_noise<-c(K_noise,z$comp$cor_noise[2])
      P_noise<-c(P_noise,z$comp$cor_noise[3])
      S_pop<-c(S_pop,z$comp$cor_pop[1])
      K_pop<-c(K_pop,z$comp$cor_pop[2])
      P_pop<-c(P_pop,z$comp$cor_pop[3])
      
      
      temp_Corl_noise<-Corbds(vi = z$last_num_keep_noise$last_num_keep_noise_c[,1],
                              vj = z$last_num_keep_noise$last_num_keep_noise_c[,2],lb = lb,ub = ub)
      Corl_noise<-c(Corl_noise,temp_Corl_noise)
      
      temp_Corl_pop<-Corbds(vi = z$last_num_keep_pop$last_num_keep_pop_c[,1],
                              vj = z$last_num_keep_pop$last_num_keep_pop_c[,2],lb = lb,ub = ub)
      
      #cat("numsim=",i,"corl_pop=",temp_Corl_pop,"\n")
      #if(is.na(temp_Corl_pop)==T){
      #  plot(z$last_num_keep_pop$last_num_keep_pop_c[,1],z$last_num_keep_pop$last_num_keep_pop_c[,2],asp=1,xlim=c(0,1),ylim=c(0,1),cex=0.5,col="red",pch=20)
      #  mtext(paste0("bad numsim = ",i,", corcoef = ",corcoef,sep=""),line=0.1)
      #  abline(a=2*lb,b=-1,col="blue")
      #  abline(a=2*ub,b=-1,col="green4")
      #  abline(a=0,b=1)
      #  rect(0,0,1,1)
      #}
     
      Corl_pop<-c(Corl_pop,temp_Corl_pop)
      
      temp_Coru_noise<-Corbds(vi = z$last_num_keep_noise$last_num_keep_noise_c[,1],
                              vj = z$last_num_keep_noise$last_num_keep_noise_c[,2],lb = 1-ub,ub = 1-lb)
      Coru_noise<-c(Coru_noise,temp_Coru_noise)
      
      temp_Coru_pop<-Corbds(vi = z$last_num_keep_pop$last_num_keep_pop_c[,1],
                            vj = z$last_num_keep_pop$last_num_keep_pop_c[,2],lb = 1-ub,ub = 1-lb)
      Coru_pop<-c(Coru_pop,temp_Coru_pop)
      
      CorlmCoru_noise<-c(CorlmCoru_noise,temp_Corl_noise-temp_Coru_noise)
      CorlmCoru_pop<-c(CorlmCoru_pop,temp_Corl_pop-temp_Coru_pop)
      
      temp_Pl_noise<-Pbds(vi = z$last_num_keep_noise$last_num_keep_noise_c[,1],
                          vj = z$last_num_keep_noise$last_num_keep_noise_c[,2],lb = lb,ub = ub)$abs_res
      Pl_noise<-c(Pl_noise,temp_Pl_noise)
      
      temp_Pl_pop<-Pbds(vi = z$last_num_keep_pop$last_num_keep_pop_c[,1],
                        vj = z$last_num_keep_pop$last_num_keep_pop_c[,2],lb = lb,ub = ub)$abs_res
      Pl_pop<-c(Pl_pop,temp_Pl_pop)
      
      temp_Pu_noise<-Pbds(vi = z$last_num_keep_noise$last_num_keep_noise_c[,1],
                          vj = z$last_num_keep_noise$last_num_keep_noise_c[,2],lb = 1-ub,ub = 1-lb)$abs_res
      Pu_noise<-c(Pu_noise,temp_Pu_noise)
      
      temp_Pu_pop<-Pbds(vi = z$last_num_keep_pop$last_num_keep_pop_c[,1],
                        vj = z$last_num_keep_pop$last_num_keep_pop_c[,2],lb = 1-ub,ub = 1-lb)$abs_res
      Pu_pop<-c(Pu_pop,temp_Pu_pop)
      
      PlmPu_noise<-c(PlmPu_noise,temp_Pl_noise-temp_Pu_noise)
      PlmPu_pop<-c(PlmPu_pop,temp_Pl_pop-temp_Pu_pop)
      
      temp_D2u_noise<-D2bds(vi = z$last_num_keep_noise$last_num_keep_noise_c[,1],
                            vj = z$last_num_keep_noise$last_num_keep_noise_c[,2],lb = 1-ub,ub = 1-lb)
      D2u_noise<-c(D2u_noise,temp_D2u_noise)
      
      temp_D2u_pop<-D2bds(vi = z$last_num_keep_pop$last_num_keep_pop_c[,1],
                          vj = z$last_num_keep_pop$last_num_keep_pop_c[,2],lb = 1-ub,ub = 1-lb)
      D2u_pop<-c(D2u_pop,temp_D2u_pop)
      
      temp_D2l_noise<-D2bds(vi = z$last_num_keep_noise$last_num_keep_noise_c[,1],
                            vj = z$last_num_keep_noise$last_num_keep_noise_c[,2],lb = lb,ub = ub)
      D2l_noise<-c(D2l_noise,temp_D2l_noise)
      
      temp_D2l_pop<-D2bds(vi = z$last_num_keep_pop$last_num_keep_pop_c[,1],
                          vj = z$last_num_keep_pop$last_num_keep_pop_c[,2],lb = lb,ub = ub)
      D2l_pop<-c(D2l_pop,temp_D2l_pop)
      
      D2umD2l_noise<-c(D2umD2l_noise,temp_D2u_noise-temp_D2l_noise)
      D2umD2l_pop<-c(D2umD2l_pop,temp_D2u_pop-temp_D2l_pop)
      
    }
    
    
    if(model=="ricker"){
      
      tempo<-paste(resloc,"hist_pop",sep="")
      if (!dir.exists(tempo)){
        dir.create(tempo)
      }
      
      K<-params[2]
      
      # Plot of histogram of pop_q (last num_keep_last number of points) from all numsim simulations
      pdf(paste0(tempo,"/",BiCopName(fcode,short=T),"_method_",method,"_corcoef_",corcoef,".pdf",sep=""),height=6,width=12)
      op<-par(mfrow=c(1,2),mar=c(5.2,4.2,1.2,1.2))
      
      skw1<-skewness(all_numsim_pop_patch1,type=2)
      hist(all_numsim_pop_patch1,breaks=1000,main=paste("patch1, skewness =",round(skw1,4),sep=""),xlab="pop. from all simulations")
      abline(v=K,col="red")
      
      skw2<-skewness(all_numsim_pop_patch2,type=2)
      hist(all_numsim_pop_patch2,breaks=1000,main=paste("patch2, skewness =",round(skw2,4),sep=""),xlab="pop. from all simulations")
      abline(v=K,col="red")
       
      par(op)
      dev.off()
      
      # Plot of time series for pop_q (last num_keep_last number of points) from each of total numsim number of simulations
      pts_patch1<-split(all_numsim_pop_patch1,as.numeric(gl(length(all_numsim_pop_patch1),num_keep_last,length(all_numsim_pop_patch1)))) 
      pts_patch2<-split(all_numsim_pop_patch1,as.numeric(gl(length(all_numsim_pop_patch1),num_keep_last,length(all_numsim_pop_patch1)))) 
      xlb<-paste("last ",num_keep_last," time series",sep="")
      tt<-c(1:num_keep_last)
      
      pdf(paste0(tempo,"/",BiCopName(fcode,short=T),"_method_",method,"_corcoef_",corcoef,"_poptimeseries_patch1.pdf",sep=""),height=25,width=25)
      op<-par(mfrow=c(5,5),mar=c(5.2,4.2,2,1.2))
      
      for(i in c(1:length(pts_patch1))){
        plot(tt,pts_patch1[[i]],xlab=xlb,ylab="pop",main=paste("numsim =",i,sep=""),type="l")
        abline(h=K,col="red")
      }
      
      par(op)
      dev.off()
      
      pdf(paste0(tempo,"/",BiCopName(fcode,short=T),"_method_",method,"_corcoef_",corcoef,"_poptimeseries_patch2.pdf",sep=""),height=25,width=25)
      op<-par(mfrow=c(5,5),mar=c(5.2,4.2,2,1.2))
      
      for(i in c(1:length(pts_patch2))){
        plot(tt,pts_patch2[[i]],xlab=xlb,ylab="pop",main=paste("numsim =",i,sep=""),type="l")
        abline(h=K,col="red")
      }
      
      par(op)
      dev.off()
      
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
    CorlmCoru_noise_mat[counter,]<-MCI(CorlmCoru_noise)
    CorlmCoru_pop_mat[counter,]<-MCI(CorlmCoru_pop)
    Pl_noise_mat[counter,]<-MCI(Pl_noise)
    Pl_pop_mat[counter,]<-MCI(Pl_pop)
    Pu_noise_mat[counter,]<-MCI(Pu_noise)
    Pu_pop_mat[counter,]<-MCI(Pu_pop)
    PlmPu_noise_mat[counter,]<-MCI(PlmPu_noise)
    PlmPu_pop_mat[counter,]<-MCI(PlmPu_pop)
    D2u_noise_mat[counter,]<-MCI(D2u_noise)
    D2u_pop_mat[counter,]<-MCI(D2u_pop)
    D2l_noise_mat[counter,]<-MCI(D2l_noise)
    D2l_pop_mat[counter,]<-MCI(D2l_pop)
    D2umD2l_noise_mat[counter,]<-MCI(D2umD2l_noise)
    D2umD2l_pop_mat[counter,]<-MCI(D2umD2l_pop)
    
    pval_S<-c(pval_S,t.test(S_noise,S_pop,alternative="two.sided",paired=T)$p.value)
    pval_K<-c(pval_K,t.test(K_noise,K_pop,alternative="two.sided",paired=T)$p.value)
    pval_P<-c(pval_P,t.test(P_noise,P_pop,alternative="two.sided",paired=T)$p.value)
    pval_Corl<-c(pval_Corl,t.test(Corl_noise,Corl_pop,alternative="two.sided",paired=T)$p.value)
    pval_Coru<-c(pval_Coru,t.test(Coru_noise,Coru_pop,alternative="two.sided",paired=T)$p.value)
    pval_CorlmCoru<-c(pval_CorlmCoru,t.test(CorlmCoru_noise,CorlmCoru_pop,alternative="two.sided",paired=T)$p.value)
    pval_Pl<-c(pval_Pl,t.test(Pl_noise,Pl_pop,alternative="two.sided",paired=T)$p.value)
    pval_Pu<-c(pval_Pu,t.test(Pu_noise,Pu_pop,alternative="two.sided",paired=T)$p.value)
    pval_PlmPu<-c(pval_PlmPu,t.test(PlmPu_noise,PlmPu_pop,alternative="two.sided",paired=T)$p.value)
    pval_D2u<-c(pval_D2u,t.test(D2u_noise,D2u_pop,alternative="two.sided",paired=T)$p.value)
    pval_D2l<-c(pval_D2l,t.test(D2l_noise,D2l_pop,alternative="two.sided",paired=T)$p.value)
    pval_D2umD2l<-c(pval_D2umD2l,t.test(D2umD2l_noise,D2umD2l_pop,alternative="two.sided",paired=T)$p.value)
    
  }
  
  return(list(resloc=resloc,
              fcode=fcode,
              method=method,
              corcoef_list=corcoef_list,
              S_noise_mat=S_noise_mat,
              K_noise_mat=K_noise_mat,
              P_noise_mat=P_noise_mat,
              S_pop_mat=S_pop_mat,
              K_pop_mat=K_pop_mat,
              P_pop_mat=P_pop_mat,
              Corl_noise_mat=Corl_noise_mat,
              Corl_pop_mat=Corl_pop_mat,
              Coru_noise_mat=Coru_noise_mat,
              Coru_pop_mat=Coru_pop_mat,
              CorlmCoru_noise_mat=CorlmCoru_noise_mat,
              CorlmCoru_pop_mat=CorlmCoru_pop_mat,
              Pl_noise_mat=Pl_noise_mat,
              Pl_pop_mat=Pl_pop_mat,
              Pu_noise_mat=Pu_noise_mat,
              Pu_pop_mat=Pu_pop_mat,
              PlmPu_noise_mat=PlmPu_noise_mat,
              PlmPu_pop_mat=PlmPu_pop_mat,
              D2u_noise_mat=D2u_noise_mat,
              D2u_pop_mat=D2u_pop_mat,
              D2l_noise_mat=D2l_noise_mat,
              D2l_pop_mat=D2l_pop_mat,
              D2umD2l_noise_mat=D2umD2l_noise_mat,
              D2umD2l_pop_mat=D2umD2l_pop_mat,
              pval_S=pval_S,
              pval_K=pval_K,
              pval_P=pval_P,
              pval_Corl=pval_Corl,
              pval_Coru=pval_Coru,
              pval_CorlmCoru=pval_CorlmCoru,
              pval_Pl=pval_Pl,
              pval_Pu=pval_Pu,
              pval_PlmPu=pval_PlmPu,
              pval_D2u=pval_D2u,
              pval_D2l=pval_D2l,
              pval_D2umD2l=pval_D2umD2l))
}
#-----------------------------------------------

