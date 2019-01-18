library(copula)
library(VineCopula)
source("CopulaFunctions_flexible.R")
source("MyBiCopGofTest.R")
#---------------------------------------------------------------------------------------
# This function GetNoise generates noises with a copula structure as you specified
# Input :
#        N : number of points drawn for a copula (C,G,J,F,SC,SG,SJ)
#        fcode : familycode of the desired copula [within c(3:6,13,14,16)]
#        corcoef : a number which may be Kendall's Tau or Spearman's Rho
#        method : a character of spearman or kendall
#        ploton : logical to genarate an optional plot
# Output :
#       list of 3 : 
#                  noise_c : a N by 2 noise matrix ,
#                  noise_q : it's qnorm transformed form 
#                  param : parameter of the copula from where that noise is generated
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
#------------------------------------------------------------------------------
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
  
  for (counter in 2:(N+1)){
    if(model=="ar1"){
      beta<-params
      alpha<-sqrt(1-beta^2)
      res[counter,]<-(beta*res[counter-1,])+(alpha*noise[counter-1,])
    }else if(model=="ricker"){
      r<-params[1]
      K<-params[2]
      res[counter,]<-res[counter-1,]*exp((r*(1-(res[counter-1,]/K)))+noise[counter-1,])
    }else{
      warning("model not specified",immediate.=T,call.=T)
    }
    
  }
  
  res2<-pnorm(res) # convert into copula space
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
  m<-mean(x)
  se<-sd(x)/sqrt(length(x))   
  return(c(m-1.96*se,m,m+1.96*se))
}
#------------------------------------------------------------------------------------------------------------
#function to get a comparison table 
# Input:
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
#s<-GetNoise(N=2000,fcode=3,corcoef=0.8,method="spearman",ploton=T)
#s2<-Simulator_Cause4copula(cons=0.5,p0=c(0,0),noise=s$noise_q)
#hist(s2$pop_q[,1],breaks=1000) # it's normal
#hist(s2$pop_c[,1],breaks=1000) # it's uniform
#points(s2$pop_c[,1],s2$pop_c[,2],col="red")
#zz<-comp(s=s,s2=s2)
#---------------------------------------------------------
#-----------------------------------------------------------------------------------------------
# This function will generate a plot for C, SC, G and SG copula (or any other as the case may be) : 
#  parameter of noise copula vs. correln coef (spearman/kendall)
# BS : number of bootstrapps used for BiCopGOFTest

Plotter_Cause4copula_GOF<-function(N,fcode,method,num_keep_last,BS,params,p0,model){
  
  corcoef_list<-seq(from=0.1,to=0.9,by=0.1)
  
  # initialize
  par_noise<-c()
  par_pop<-c()
  se_par_pop<-c()
  pval_CvM<-c()
  pval_KS<-c()
  BS_success<-c()
  
 
  for(corcoef in corcoef_list){
    
    #cat("corcoef=",corcoef,"\n")
    
    # generate noise copula 
    s<-GetNoise(N=N,fcode=fcode,corcoef=corcoef,method=method,ploton=F)
    # generate pop_copula
    s2<-Simulator_Cause4copula(params=params,p0=cp0,noise=s$noise_q,model=model)
    
    # compare btw parameters of s$noise_c and s2$pop_c
    
    u1<-tail(s2$pop_c[,1],num_keep_last)  # take last num_keep_last # of rows of s2$pop_c matrix
    u2<-tail(s2$pop_c[,2],num_keep_last)
    
    z<-BiCopEst(u1,u2,family=fcode,se=T) # apply BiCopEst to get par_pop and se_par_pop for sample estimation
    zf<-MyBiCopGofTest(u1,u2,family=fcode,method = "kendall",B=BS) # to get p-values(CvM, KS) of GOF test
    
    par_noise<-c(par_noise,s$param)
    par_pop<-c(par_pop,z$par)
    se_par_pop<-c(se_par_pop,z$se)
    pval_CvM<-c(pval_CvM,zf$p.value.CvM)
    pval_KS<-c(pval_KS,zf$p.value.KS)
    BS_success<-c(BS_success,zf$B_success)
    
  }
  
 
  BS_success_percentage<-((min(BS_success))/BS)*100
  
  if(method=="spearman"){
    xlabel<-"Spearman's Rho"
  }else if(method=="kendall"){
    xlabel<-"Kendall's Tau"
  }else{
    warning("specify method",immediate.=T,call.=T)
  }
  
    ## add extra space to right margin of plot within frame
    op<-par(mar=c(3.5, 4, 4, 6) + 0.1)
    se_par_pop_lim<-max(se_par_pop,na.rm=T) # to remove NA from ylim
    ## Plot first set of data and draw its axis
    plot(corcoef_list, par_pop, pch=6, axes=FALSE, 
         #ylim=c(ceiling(min(0,par_pop-1.96*se_par_pop_lim)),
         ylim=c(0,ceiling(max(par_noise,par_pop+1.96*se_par_pop_lim))),xlim=c(0,1),
         xlab="", ylab="", 
         col="blue")#,main=BiCopName(family = fcode,short=F))
    segments(corcoef_list,par_pop-1.96*se_par_pop,corcoef_list,par_pop+1.96*se_par_pop,col='blue')
    bar_len<-0.02
    segments(corcoef_list-bar_len,par_pop-1.96*se_par_pop,corcoef_list+bar_len,par_pop-1.96*se_par_pop,col='blue')
    segments(corcoef_list-bar_len,par_pop+1.96*se_par_pop,corcoef_list+bar_len,par_pop+1.96*se_par_pop,col='blue')
    #arrows(corcoef_list,par_pop-1.96*se_par_pop,corcoef_list,par_pop+1.96*se_par_pop,length=0.03, angle=90, code=3, col='blue')
    points(corcoef_list,par_noise,pch=2,col="red")
    axis(2, col="black",las=1)  ## las=1 makes horizontal labels
    mtext(expression(paste("Parameter (",theta,")")), side = 2, line = 2.5)
    mtext(paste0(BiCopName(family = fcode,short=T)," , min_BS_success = ",BS_success_percentage,"%"),side=3,line=0.3,cex=0.8)
    box()
    
    ## Allow a second plot on the same graph
    par(new=TRUE)
    
    ## Plot the second plot and put axis scale on right
    plot(corcoef_list, pval_CvM, pch=16,  xlab="", ylab="", ylim=c(0,1), xlim=c(0,1),
         axes=FALSE, type="p", col="magenta")
    points(corcoef_list, pval_KS, pch=16, col="green2")
    lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
    ## a little farther out (line=3) to make room for labels
    mtext("p values",side=4,col="purple",line=3) 
    axis(4, ylim=c(0,1), col="purple",col.axis="purple",las=1)
    
    ## Draw the x axis
    axis(1,corcoef_list)
    mtext(xlabel,side=1,col="black",line=2.5)  
    par(op)

}

#--------------------------------------------------------------------------------------
# This function will give you a 3 by 3 multipanel plot for a specified family of copula 
# 1st row : [1] Spearman correlation btw noise_c & pop_c against spearman's rho/ kendall's tau of noise
#           [2] Kendall correlation btw noise_c & pop_c against spearman's rho/ kendall's tau of noise
#           [3] Pearson correlation btw noise_q & pop_q against spearman's rho/ kendall's tau of noise
# 2nd and 3rd row :
#  Our non-parametric stats @ extreme ends (Corl,Coru,Pl,Pu,D2u,D2l) 
#                                 against spearman's rho/ kendall's tau of noise
#-------------------------------------------------------------------------------------
#
# Input :
#       numsim : a number over which desired stat (Spearman, Kendall, Pearson) called for (default:50)
#       fcode : family of copula from where noise is genarated initially [within c(3:6,13,14,16)]
#       method : a character : either "spearman" or "kendall"
#       lb : lower bound for Non-parametric stat function (Default=0)
#       ub : upper bound for Non-parametric stat function (Default =0.1)
#       num_keep_last : number of rows you want to keep from bottom 
#                     for each of s$noise_c, s$noise_q and s2$pop_c, s2$pop_q matrix
#       resloc : folder location where the plots should be saved
#       params,p0,model : these are inputs for Simulator_Cause4copula function
#-----------------------------------------------------------------------------------------------
Plotter_Cause4copula_stat<-function(N,numsim=50,fcode,method,lb=0,ub=0.1,num_keep_last,resloc,params,p0,model){
  
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
    
    for(i in 1:numsim){
      #cat("i=",i,"\n")
      s<-GetNoise(N=N,fcode=fcode,corcoef=corcoef,method=method,ploton=F)
      s2<-Simulator_Cause4copula(params=params,p0=p0,noise=s$noise_q,model=model)
      z<-comp(s=s,s2=s2,num_keep_last = num_keep_last)
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
  
  if(method=="spearman"){
    xlabel<-"Spearman's Rho"
  }else if(method=="kendall"){
    xlabel<-"Kendall's Tau"
  }else{
    warning("specify method",immediate.=T,call.=T)
  }
 
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Spearman_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,S_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Spearman",xlim=c(0,1),ylim=c(0,1),cex.lab=2,cex.axis=1.5)
  segments(corcoef_list,S_noise_mat[,1],corcoef_list,S_noise_mat[,3],col='red')
  bar_len<-0.02
  segments(corcoef_list-bar_len,S_noise_mat[,1],corcoef_list+bar_len,S_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,S_noise_mat[,3],corcoef_list+bar_len,S_noise_mat[,3],col='red')
  #arrows(corcoef_list,S_noise_mat[,1],corcoef_list,S_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,S_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,S_pop_mat[,1],corcoef_list,S_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,S_pop_mat[,1],corcoef_list+bar_len,S_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,S_pop_mat[,3],corcoef_list+bar_len,S_pop_mat[,3],col='blue')
  #arrows(corcoef_list,S_pop_mat[,1],corcoef_list,S_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  points(corcoef_list,pval_S,col="purple")
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Kendall_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,K_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Kendall",xlim=c(0,1),ylim=c(0,1),cex.lab=2,cex.axis=1.5)
  segments(corcoef_list,K_noise_mat[,1],corcoef_list,K_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,K_noise_mat[,1],corcoef_list+bar_len,K_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,K_noise_mat[,3],corcoef_list+bar_len,K_noise_mat[,3],col='red')
  #arrows(corcoef_list,K_noise_mat[,1],corcoef_list,K_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,K_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,K_pop_mat[,1],corcoef_list,K_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,K_pop_mat[,1],corcoef_list+bar_len,K_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,K_pop_mat[,3],corcoef_list+bar_len,K_pop_mat[,3],col='blue')
  #arrows(corcoef_list,K_pop_mat[,1],corcoef_list,K_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  points(corcoef_list,pval_K,col="purple")
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pearson_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,P_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab="Pearson",xlim=c(0,1),ylim=c(0,1),cex.lab=2,cex.axis=1.5)
  segments(corcoef_list,P_noise_mat[,1],corcoef_list,P_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,P_noise_mat[,1],corcoef_list+bar_len,P_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,P_noise_mat[,3],corcoef_list+bar_len,P_noise_mat[,3],col='red')
  #arrows(corcoef_list,P_noise_mat[,1],corcoef_list,P_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,P_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,P_pop_mat[,1],corcoef_list,P_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,P_pop_mat[,1],corcoef_list+bar_len,P_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,P_pop_mat[,3],corcoef_list+bar_len,P_pop_mat[,3],col='blue')
  #arrows(corcoef_list,P_pop_mat[,1],corcoef_list,P_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  points(corcoef_list,pval_P,col="purple")
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Corl_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Corl_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("Cor"["l"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.1+max(Corl_noise_mat[,2],Corl_pop_mat[,2])))
  segments(corcoef_list,Corl_noise_mat[,1],corcoef_list,Corl_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,Corl_noise_mat[,1],corcoef_list+bar_len,Corl_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,Corl_noise_mat[,3],corcoef_list+bar_len,Corl_noise_mat[,3],col='red')
  #arrows(corcoef_list,Corl_noise_mat[,1],corcoef_list,Corl_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Corl_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,Corl_pop_mat[,1],corcoef_list,Corl_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,Corl_pop_mat[,1],corcoef_list+bar_len,Corl_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,Corl_pop_mat[,3],corcoef_list+bar_len,Corl_pop_mat[,3],col='blue')
  #arrows(corcoef_list,Corl_pop_mat[,1],corcoef_list,Corl_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Corl,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Coru_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Coru_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("Cor"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.1+max(Coru_noise_mat[,2],Coru_pop_mat[,2])))
  segments(corcoef_list,Coru_noise_mat[,1],corcoef_list,Coru_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,Coru_noise_mat[,1],corcoef_list+bar_len,Coru_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,Coru_noise_mat[,3],corcoef_list+bar_len,Coru_noise_mat[,3],col='red')
  #arrows(corcoef_list,Coru_noise_mat[,1],corcoef_list,Coru_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Coru_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,Coru_pop_mat[,1],corcoef_list,Coru_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,Coru_pop_mat[,1],corcoef_list+bar_len,Coru_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,Coru_pop_mat[,3],corcoef_list+bar_len,Coru_pop_mat[,3],col='blue')
  #arrows(corcoef_list,Coru_pop_mat[,1],corcoef_list,Coru_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Coru,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Corl-Coru_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,CorlmCoru_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("Cor"["l"]-"Cor"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(-0.2,0.2))
  segments(corcoef_list,CorlmCoru_noise_mat[,1],corcoef_list,CorlmCoru_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,CorlmCoru_noise_mat[,1],corcoef_list+bar_len,CorlmCoru_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,CorlmCoru_noise_mat[,3],corcoef_list+bar_len,CorlmCoru_noise_mat[,3],col='red')
  #arrows(corcoef_list,CorlmCoru_noise_mat[,1],corcoef_list,CorlmCoru_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,CorlmCoru_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,CorlmCoru_pop_mat[,1],corcoef_list,CorlmCoru_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,CorlmCoru_pop_mat[,1],corcoef_list+bar_len,CorlmCoru_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,CorlmCoru_pop_mat[,3],corcoef_list+bar_len,CorlmCoru_pop_mat[,3],col='blue')
  #arrows(corcoef_list,CorlmCoru_pop_mat[,1],corcoef_list,CorlmCoru_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  lines(range(0,1),c(0,0),type='l',lty='dotted',col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_CorlmCoru,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_scatter_CorlmCoru_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(4.5,6,2,1), mgp=c(3,0.5,0))
  ylim1<-range(CorlmCoru_pop_mat[,1],CorlmCoru_pop_mat[,3])
  xlim1<-range(CorlmCoru_noise_mat[,1],CorlmCoru_noise_mat[,3])
  b_len<-diff(xlim1)/50
  plot(CorlmCoru_noise_mat[,2],CorlmCoru_pop_mat[,2],cex=2,col="black",
       xlab=expression("Cor"["l,noise"]-"Cor"["u,noise"]),
       ylab=expression("Cor"["l,pop"]-"Cor"["u,pop"]),
       cex.lab=2,cex.axis=1.5,
       xlim=xlim1,ylim=ylim1)
  segments(CorlmCoru_noise_mat[,2],CorlmCoru_pop_mat[,1],CorlmCoru_noise_mat[,2],CorlmCoru_pop_mat[,3],col='blue')
  segments(CorlmCoru_noise_mat[,2]-b_len,CorlmCoru_pop_mat[,1],CorlmCoru_noise_mat[,2]+b_len,CorlmCoru_pop_mat[,1],col='blue')
  segments(CorlmCoru_noise_mat[,2]-b_len,CorlmCoru_pop_mat[,3],CorlmCoru_noise_mat[,2]+b_len,CorlmCoru_pop_mat[,3],col='blue')
  
  segments(CorlmCoru_noise_mat[,1],CorlmCoru_pop_mat[,2],CorlmCoru_noise_mat[,3],CorlmCoru_pop_mat[,2],col='red')
  segments(CorlmCoru_noise_mat[,1],CorlmCoru_pop_mat[,2]-b_len,CorlmCoru_noise_mat[,1],CorlmCoru_pop_mat[,2]+b_len,col='red')
  segments(CorlmCoru_noise_mat[,3],CorlmCoru_pop_mat[,2]-b_len,CorlmCoru_noise_mat[,3],CorlmCoru_pop_mat[,2]+b_len,col='red')
  
  dat<-data.frame(x1=CorlmCoru_noise_mat[,2],y1=CorlmCoru_pop_mat[,2])
  mylm<-lm(y1~x1,data=dat)
  abline(mylm,col="black")
  lines(x=xlim1,y=ylim1,col="green")
  c<-cor.test(dat$x1,dat$y1,method = "pearson",alternative = "t")
  mtext(paste0("Pearson correlation = ",round(unname(c$estimate),3),", p = ",round(c$p.value,4),sep=""),cex=1.5,line=0.1)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pl_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Pl_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("P"["l"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.02+max(Pl_noise_mat[,2],Pl_pop_mat[,2])))
  segments(corcoef_list,Pl_noise_mat[,1],corcoef_list,Pl_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,Pl_noise_mat[,1],corcoef_list+bar_len,Pl_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,Pl_noise_mat[,3],corcoef_list+bar_len,Pl_noise_mat[,3],col='red')
  #arrows(corcoef_list,Pl_noise_mat[,1],corcoef_list,Pl_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Pl_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,Pl_pop_mat[,1],corcoef_list,Pl_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,Pl_pop_mat[,1],corcoef_list+bar_len,Pl_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,Pl_pop_mat[,3],corcoef_list+bar_len,Pl_pop_mat[,3],col='blue')
  #arrows(corcoef_list,Pl_pop_mat[,1],corcoef_list,Pl_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Pl,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  #op<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 1, 0, 0), new = TRUE)
  #plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #legend("top", c("noise copula","population copula", "p-value(paired t-test)"), col = c("red", "blue", "purple"),
  #       cex = 0.8, pch = c(16, 16, 1), xpd = TRUE, horiz = TRUE, inset = c(0,0),
  #       bty = "n") 
  #par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pu_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Pu_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("P"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.02+max(Pu_noise_mat[,2],Pu_pop_mat[,2])))
  segments(corcoef_list,Pu_noise_mat[,1],corcoef_list,Pu_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,Pu_noise_mat[,1],corcoef_list+bar_len,Pu_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,Pu_noise_mat[,3],corcoef_list+bar_len,Pu_noise_mat[,3],col='red')
  #arrows(corcoef_list,Pu_noise_mat[,1],corcoef_list,Pu_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,Pu_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,Pu_pop_mat[,1],corcoef_list,Pu_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,Pu_pop_mat[,1],corcoef_list+bar_len,Pu_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,Pu_pop_mat[,3],corcoef_list+bar_len,Pu_pop_mat[,3],col='blue')
  #arrows(corcoef_list,Pu_pop_mat[,1],corcoef_list,Pu_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_Pu,col="purple",type="p",axes = FALSE,
       bty = "n", xlab = "", ylab = "",ylim=c(0,1),xlim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pl-Pu_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,PlmPu_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("P"["l"]-"P"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(-0.1,0.1))
  segments(corcoef_list,PlmPu_noise_mat[,1],corcoef_list,PlmPu_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,PlmPu_noise_mat[,1],corcoef_list+bar_len,PlmPu_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,PlmPu_noise_mat[,3],corcoef_list+bar_len,PlmPu_noise_mat[,3],col='red')
  #arrows(corcoef_list,PlmPu_noise_mat[,1],corcoef_list,PlmPu_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,PlmPu_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,PlmPu_pop_mat[,1],corcoef_list,PlmPu_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,PlmPu_pop_mat[,1],corcoef_list+bar_len,PlmPu_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,PlmPu_pop_mat[,3],corcoef_list+bar_len,PlmPu_pop_mat[,3],col='blue')
  #arrows(corcoef_list,PlmPu_pop_mat[,1],corcoef_list,PlmPu_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  lines(range(0,1),c(0,0),type='l',lty='dotted',col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_PlmPu,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_scatter_PlmPu_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(4.5,6,2,1), mgp=c(3,0.5,0))
  ylim1<-range(PlmPu_pop_mat[,1],PlmPu_pop_mat[,3])
  xlim1<-range(PlmPu_noise_mat[,1],PlmPu_noise_mat[,3])
  
  b_len<-diff(xlim1)/50
  
  plot(PlmPu_noise_mat[,2],PlmPu_pop_mat[,2],cex=2,col="black",
       xlab=expression("P"["l,noise"]-"P"["u,noise"]),
       ylab=expression("P"["l,pop"]-"P"["u,pop"]),
       cex.lab=2,cex.axis=1.5,
       xlim=xlim1,ylim=ylim1)
  
  segments(PlmPu_noise_mat[,2],PlmPu_pop_mat[,1],PlmPu_noise_mat[,2],PlmPu_pop_mat[,3],col='blue')
  segments(PlmPu_noise_mat[,2]-b_len,PlmPu_pop_mat[,1],PlmPu_noise_mat[,2]+b_len,PlmPu_pop_mat[,1],col='blue')
  segments(PlmPu_noise_mat[,2]-b_len,PlmPu_pop_mat[,3],PlmPu_noise_mat[,2]+b_len,PlmPu_pop_mat[,3],col='blue')
  
  segments(PlmPu_noise_mat[,1],PlmPu_pop_mat[,2],PlmPu_noise_mat[,3],PlmPu_pop_mat[,2],col='red')
  segments(PlmPu_noise_mat[,1],PlmPu_pop_mat[,2]-b_len,PlmPu_noise_mat[,1],PlmPu_pop_mat[,2]+b_len,col='red')
  segments(PlmPu_noise_mat[,3],PlmPu_pop_mat[,2]-b_len,PlmPu_noise_mat[,3],PlmPu_pop_mat[,2]+b_len,col='red')
  
  dat<-data.frame(x1=PlmPu_noise_mat[,2],y1=PlmPu_pop_mat[,2])
  mylm<-lm(y1~x1,data=dat)
  abline(mylm,col="black")
  lines(x=xlim1,y=ylim1,col="green")
  c<-cor.test(dat$x1,dat$y1,method = "pearson",alternative = "t")
  mtext(paste0("Pearson correlation = ",round(unname(c$estimate),3),", p = ",round(c$p.value,4),sep=""),cex=1.5,line=0.1)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_D2u_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,D2u_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("D"[u]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.002+max(D2u_noise_mat[,2],D2u_pop_mat[,2])))
  segments(corcoef_list,D2u_noise_mat[,1],corcoef_list,D2u_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,D2u_noise_mat[,1],corcoef_list+bar_len,D2u_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,D2u_noise_mat[,3],corcoef_list+bar_len,D2u_noise_mat[,3],col='red')
  #arrows(corcoef_list,D2u_noise_mat[,1],corcoef_list,D2u_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,D2u_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,D2u_pop_mat[,1],corcoef_list,D2u_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,D2u_pop_mat[,1],corcoef_list+bar_len,D2u_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,D2u_pop_mat[,3],corcoef_list+bar_len,D2u_pop_mat[,3],col='blue')
  #arrows(corcoef_list,D2u_pop_mat[,1],corcoef_list,D2u_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_D2u,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_D2l_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,D2l_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("D"[l]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.002+max(D2l_noise_mat[,2],D2l_pop_mat[,2])))
  segments(corcoef_list,D2l_noise_mat[,1],corcoef_list,D2l_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,D2l_noise_mat[,1],corcoef_list+bar_len,D2l_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,D2l_noise_mat[,3],corcoef_list+bar_len,D2l_noise_mat[,3],col='red')
  #arrows(corcoef_list,D2l_noise_mat[,1],corcoef_list,D2l_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,D2l_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,D2l_pop_mat[,1],corcoef_list,D2l_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,D2l_pop_mat[,1],corcoef_list+bar_len,D2l_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,D2l_pop_mat[,3],corcoef_list+bar_len,D2l_pop_mat[,3],col='blue')
  #arrows(corcoef_list,D2l_pop_mat[,1],corcoef_list,D2l_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  par(new = TRUE)
  plot(corcoef_list,pval_D2l,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_D2u-D2l_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,D2umD2l_noise_mat[,2],cex=0.5,col="red",xlab=xlabel,ylab=expression("D"["u"]^2-"D"["l"]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(-0.01,0.01)) 
       #ylim=c(0,0.002+max(D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,2])))
  segments(corcoef_list,D2umD2l_noise_mat[,1],corcoef_list,D2umD2l_noise_mat[,3],col='red')
  segments(corcoef_list-bar_len,D2umD2l_noise_mat[,1],corcoef_list+bar_len,D2umD2l_noise_mat[,1],col='red')
  segments(corcoef_list-bar_len,D2umD2l_noise_mat[,3],corcoef_list+bar_len,D2umD2l_noise_mat[,3],col='red')
  #arrows(corcoef_list,D2umD2l_noise_mat[,1],corcoef_list,D2umD2l_noise_mat[,3],length=0.03, angle=90, code=3, col='red')
  points(corcoef_list,D2umD2l_pop_mat[,2],cex=0.5,col="blue")
  segments(corcoef_list,D2umD2l_pop_mat[,1],corcoef_list,D2umD2l_pop_mat[,3],col='blue')
  segments(corcoef_list-bar_len,D2umD2l_pop_mat[,1],corcoef_list+bar_len,D2umD2l_pop_mat[,1],col='blue')
  segments(corcoef_list-bar_len,D2umD2l_pop_mat[,3],corcoef_list+bar_len,D2umD2l_pop_mat[,3],col='blue')
  #arrows(corcoef_list,D2umD2l_pop_mat[,1],corcoef_list,D2umD2l_pop_mat[,3],length=0.03, angle=90, code=3, col='blue')
  lines(range(0,1),c(0,0),type='l',lty='dotted',col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_D2umD2l,col="purple",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1))
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='purple')
  axis(side=4,col='purple',col.axis="purple",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='purple',cex=2)
  par(op)
  dev.off()
  
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_scatter_D2umD2l_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(4.5,6,2,1), mgp=c(3,0.5,0))
  ylim1<-range(D2umD2l_pop_mat[,1],D2umD2l_pop_mat[,3])
  xlim1<-range(D2umD2l_noise_mat[,1],D2umD2l_noise_mat[,3])
  
  plot(D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,2],cex=2,col="black",
       xlab=expression("D"["u,noise"]^2-"D"["l,noise"]^2),
       ylab=expression("D"["u,pop"]^2-"D"["l,pop"]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=xlim1,ylim=ylim1)
  
  b_len<-diff(xlim1)/50
  segments(D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,1],D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,3],col='blue')
  segments(D2umD2l_noise_mat[,2]-b_len,D2umD2l_pop_mat[,1],D2umD2l_noise_mat[,2]+b_len,D2umD2l_pop_mat[,1],col='blue')
  segments(D2umD2l_noise_mat[,2]-b_len,D2umD2l_pop_mat[,3],D2umD2l_noise_mat[,2]+b_len,D2umD2l_pop_mat[,3],col='blue')
  
  segments(D2umD2l_noise_mat[,1],D2umD2l_pop_mat[,2],D2umD2l_noise_mat[,3],D2umD2l_pop_mat[,2],col='red')
  segments(D2umD2l_noise_mat[,1],D2umD2l_pop_mat[,2]-b_len,D2umD2l_noise_mat[,1],D2umD2l_pop_mat[,2]+b_len,col='red')
  segments(D2umD2l_noise_mat[,3],D2umD2l_pop_mat[,2]-b_len,D2umD2l_noise_mat[,3],D2umD2l_pop_mat[,2]+b_len,col='red')
  
  dat<-data.frame(x1=D2umD2l_noise_mat[,2],y1=D2umD2l_pop_mat[,2])
  mylm<-lm(y1~x1,data=dat)
  abline(mylm,col="black")
  lines(x=xlim1,y=ylim1,col="green")
  c<-cor.test(dat$x1,dat$y1,method = "pearson",alternative = "t")
  mtext(paste0("Pearson correlation = ",round(unname(c$estimate),3),", p = ",round(c$p.value,4),sep=""),cex=1.5,line=0.1)
  par(op)
  dev.off()
  
  #op2<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  #plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #legend("top", c("noise","population", "p-value(paired t-test)"), col = c("red", "blue", "purple"),
  #       cex = 0.8, pch = c(16, 16, 1), xpd = TRUE, horiz = T, inset = c(0,0), 
  #       bty = "n") 
  #par(op2)
  
  pdf(paste0(resloc,"common_legend_cause4copula_stat.pdf",sep=""),height=1,width=15)
  op<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("center", c("noise copula","population copula", "p-value(paired t-test)"), col = c("red", "blue", "purple"),
                cex = 2.5, pch = c(16, 16, 1), xpd = TRUE, horiz = T, inset = c(0,0), 
                bty = "n")
  par(op)
  dev.off()
}
#-----------------------------------------------

