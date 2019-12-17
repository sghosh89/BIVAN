#Code for comprehensively testing the stats in CopulaFunctions.R

library(copula)
library(parallel)
source("CopulaFunctions.R")
#***
#calculations to test the stats
#***

#This function is the worker for testing the stats. It generates numpts
#from cop (a copula) and then calls all the stats on those data, then repeats
#numsims times and returns all the results. The argument rank is for whether 
#you perform ranking after simulating from the copula or not.
#
#Args
#cop      A copula (2 dimensional)
#numpts   Numbers of points to pull from the copula each simulation
#numsims  The number of simulations to do and apply the statistics to
#rank     T or F, apply the ranking to sims or not?
#
#Output
#A 6 by numsims matrix with all the values of all the stats on each sim
#
worker<-function(cop,numpts,numsims,rank){
  
  res<-matrix(NA,6,numsims)
  rownames(res)<-c("Corl","Coru","Pl","Pu","D2l","D2u") 
  
  for (simcount in 1:numsims){
    #simulate points from the copula
    v<-rCopula(numpts,cop)
    vi<-v[,1]
    vj<-v[,2]
    
    #if rank==T, rank the simulated points
    if (rank==T){
      vi<-(rank(vi)/(length(vi)+1))   # gives the same value if you write as : vi<-pobs(vi)
      vj<-(rank(vj)/(length(vj)+1))   # gives the same value if you write as : vj<-pobs(vj)
    }
    
    #call each of the stats on the points and save in res
    res[1:2,simcount]<-CorlCoru(vi,vj)
    temp1<-PlPu(vi,vj)
    res[3,simcount]<-temp1[[2]]
    res[4,simcount]<-temp1[[3]]
    res[5:6,simcount]<-D2lD2u(vi,vj)
  }
  
  return(res)
}

#----------------------------------------------------------------------------------------------
#set up list of copulas on which to call worker using mclapply, start by finding 
#the copula parameters
#
# This function takes:
# Input : spearvals : a sequence of values 0 to 1 which indicates spearman correlation
# Output : a list of 2 things:
#            1. coplist : a long list (length=4*length(spearvals)) of 4 types of copula : clayton, normal, frank, survival clayton
#            2. paramlist : a long list (length=length(coplist)) of parameters of above 4 types of copula, this paramlist
#                           indicates the parameters for which each copula should have that given specified spearman correlation.
#
coplist_with_params<-function(spearvals){
  
  taCcop<-claytonCopula(3,2)  # initialize by any old Clayton and Frank and normal copulas
  taNcop<-normalCopula(.2,2,"un")
  taFcop<-frankCopula(3,2)
  
  Cparams<-NA*numeric(length(spearvals))
  Nparams<-NA*numeric(length(spearvals))
  Fparams<-NA*numeric(length(spearvals))
  
  for (counter in 1:length(spearvals))
  {
    Cparams[counter]<-iRho(taCcop,spearvals[counter])
    Nparams[counter]<-iRho(taNcop,spearvals[counter])
    Fparams[counter]<-iRho(taFcop,spearvals[counter])
  }
  
  Cparams[spearvals==0]<-0
  Fparams[spearvals==0]<-0
  Cfparams<-Cparams #for the survival (flipped) Clayton, same as the Clayton
  
  #now the copulas themselves
  coplist<-list()

  #Clayton  
  for (counter in 1:length(Cparams)){ 
    coplist<-c(coplist,claytonCopula(Cparams[counter],2))
  } 
  
  #Gaussian (normal) copula
  for (counter in 1:length(Nparams)){
    coplist<-c(coplist,normalCopula(Nparams[counter],2,"un"))
  }
  
  #Frank
  for (counter in 1:length(Fparams)){
    coplist<-c(coplist,frankCopula(Fparams[counter],2))
  }
  
  #flipped/survival Claytons
  for (counter in 1:length(Cfparams)){
    coplist<-c(coplist,rotCopula(claytonCopula(Cfparams[counter],2)))
  }
  
  paramlist<- list(Cparams=Cparams,
                   Nparams=Nparams,
                   Fparams=Fparams,
                   Cfparams=Cfparams)
  
  return(list(coplist=coplist, 
              paramlist=paramlist))
}

#--------------------------------------------------------------------------------------
# check the function
#------------------------
#callfn<-coplist_with_params(spearvals=seq(from=0,to=0.9,by=0.1))
#coplist<-callfn$coplist
#ncores<-3
#reslist<-mclapply(X=coplist,FUN=worker,numpts=35,numsims=100,rank=F,mc.cores=ncores)
#----------------------------------

# This function merges stat for (Corl-Coru), .....etc for all the simulations
# Input:
#      1. reslist : A list of 4*length(spearvals) matrices and each is a 6 by numsims matrix, 6 rows are : Corl, Coru, Pl, Pu, D2l, D2u
#      2. numsims : The number of simulations to do and apply the statistics to
#      3. j1, j2 :  1 to 6, index of those 6 stats
#
stat_list<-function(reslist,numsims,j1,j2){ 
  stat_sms<-matrix(NA,nrow=length(reslist),ncol=numsims)
  for (i in 1:length(reslist)){
    stat_sms[i,]<-(reslist[[i]][j1,]-reslist[[i]][j2,])
  } 
  return(stat_sms)
}

# This function gives mean,lowCI,upCI of a vector
#
mCI<-function(x){
  m<-mean(x)
  se<-sd(x)/sqrt(length(x))   # denominator should be sqrt(length(x)-1) ?????
  return(c(m-1.96*se,m,m+1.96*se))
}

#--------------------------------------------------------------------------------------------------
# This function gives a list of 3 matrices for (Corl-Coru), (Pl-Pu), (D2u-D2l) stat:
#                         each matrix has 3 columns [lowCI, mean, upCI] and nrows=4*length(spearvals)
#                                                   4 : type of copula(C,N,F,SC)
#                                                   length of spearvals =10 (here)
#------------------------------------------------------------------------------------------------------------------------
#
M1mM2<-function(reslist,numsims){
  
  
  stat_CorlmCoru<-stat_list(reslist,numsims,1,2)
  stat_PlmPu<-stat_list(reslist,numsims,3,4)
  stat_D2umD2l<-stat_list(reslist,numsims,6,5)
  
  CorlmCoru<-matrix(NA,nrow=dim(stat_CorlmCoru)[1],ncol=3)
  colnames(CorlmCoru)<-c('lowCI','mean','upCI')
  PlmPu<-matrix(NA,nrow=dim(stat_CorlmCoru)[1],ncol=3)
  colnames(PlmPu)<-colnames(CorlmCoru)
  D2umD2l<-matrix(NA,nrow=dim(stat_CorlmCoru)[1],ncol=3)
  colnames(D2umD2l)<-colnames(CorlmCoru)
  
  for(i in 1:dim(stat_CorlmCoru)[1]){ 
    
    temp0<-mCI(stat_CorlmCoru[i,])
    CorlmCoru[i,1]<-temp0[1]
    CorlmCoru[i,2]<-temp0[2]
    CorlmCoru[i,3]<-temp0[3]
    
    temp0<-mCI(stat_PlmPu[i,])
    PlmPu[i,1]<-temp0[1]
    PlmPu[i,2]<-temp0[2]
    PlmPu[i,3]<-temp0[3]
    
    temp0<-mCI(stat_D2umD2l[i,])
    D2umD2l[i,1]<-temp0[1]
    D2umD2l[i,2]<-temp0[2]
    D2umD2l[i,3]<-temp0[3]
    
  } 
  
  return(list(CorlmCoru=CorlmCoru,
              PlmPu=PlmPu,
              D2umD2l=D2umD2l))
}

#-----------------------------------------------------------------------------------------
# This is a t-test function which returns a p-value
#
ttest<-function(reslist,numsims,j1,j2,is,ie,ispearval){ 
  stat_M1mM2<-stat_list(reslist,numsims,j1,j2)
  x<-stat_M1mM2[is:ie,][ispearval,]
  tr<-t.test(x) # x is a vector
  t<-tr$statistic
  p<-tr$p.value
  return(p)
}

#-------------------------------------------------------------------------------------------------------
# A plotter function that takes res_pt35_rankF, etc. and makes the appropriate plot
#
plotter_stat_testing<-function(reslist,filename,xaxparams,resultsloc,spearvals){
  
  numsims<-dim(reslist[[1]])[2]
  
  a<-M1mM2(reslist,numsims)
  b<-xaxparams
  
  pdf(paste(resultsloc,filename,".pdf",sep=""),width=12, height=8)
  #op<-par(mfcol=c(3,4),mar=c(3,3.5,3,3.5), mgp=c(1.5,0.5,0))
  op<-par(mfcol=c(3,4),mar=c(4,4.5,3,4.5), mgp=c(2,0.5,0))
  
  t1<-c(1,length(b[[1]])+1,length(b[[1]])+length(b[[2]])+1,length(b[[1]])+length(b[[2]])+length(b[[3]])+1)
  t2<-c(length(b[[1]]),length(b[[1]])+length(b[[2]]),length(b[[1]])+length(b[[2]])+length(b[[3]]),length(b[[1]])+length(b[[2]])+length(b[[3]])+length(b[[4]]))
  
  for(i in c(1,2,3,4)){
    #if(i==1 | i==4){
      #is<-t1[i]
      #ie<-t2[i]
      
      #t_CorlmCoru<-ttest(reslist,numsims,1,2,is,ie,1)
     # plot(spearvals,a$CorlmCoru[t1[i]:t2[i],][,2],type='o',col='black',ylim=c(-0.25,0.25),xlab='Spearman',ylab=expression(cor[l] - cor[u]))
    #  arrows(spearvals, a$CorlmCoru[t1[i]:t2[i],][,1],spearvals,a$CorlmCoru[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='black')
     # lines(range(spearvals),c(0,0),type='l',lty='dashed',col='black')  
     # mtext(paste0("p = ",round(t_CorlmCoru,3)),side = 3, line=-1.5, adj=0.5, col='black')
      
     # t_PlmPu<-ttest(reslist,numsims,3,4,is,ie,1)
     # plot(spearvals,a$PlmPu[t1[i]:t2[i],][,2],type='o',col='black',ylim=c(-0.25,0.25),xlab='Spearman',ylab=expression(P[l] - P[u]))  
     # arrows(spearvals,a$PlmPu[t1[i]:t2[i],][,1],spearvals,a$PlmPu[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='black')
     # lines(range(spearvals),c(0,0),type='l',lty='dashed',col='black')  
     # mtext(paste0("p = ",round(t_PlmPu,3)),side = 3, line=-1.5, adj=0.5, col='black')
      
     # t_D2umD2l<-ttest(reslist,numsims,6,5,is,ie,1)
     # plot(spearvals,a$D2umD2l[t1[i]:t2[i],][,2],type='o',col='black',ylim=c(-0.05,0.05),xlab='Spearman',ylab=expression(D[u]^2 - D[l]^2)) 
     # arrows(spearvals,a$D2umD2l[t1[i]:t2[i],][,1],spearvals,a$D2umD2l[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='black')
     # lines(range(spearvals),c(0,0),type='l',lty='dashed',col='black')
     # mtext(paste0("p = ",round(t_D2umD2l,3)),side = 3, line=-1.5, adj=0.5, col='black')
      
    #}else{
      
      is<-t1[i]
      ie<-t2[i]
      
      t_CorlmCoru<-c()
      t_PlmPu<-c()
      t_D2umD2l<-c()
      
      for(ispearval in 1:length(spearvals)){
        x1<-ttest(reslist,numsims,j1=1,j2=2,is,ie,ispearval)
        x2<-ttest(reslist,numsims,j1=3,j2=4,is,ie,ispearval)
        x3<-ttest(reslist,numsims,j1=6,j2=5,is,ie,ispearval)
        
        t_CorlmCoru<-c(t_CorlmCoru,x1)
        t_PlmPu<-c(t_PlmPu,x2)
        t_D2umD2l<-c(t_D2umD2l,x3)
      }
      
      plot(spearvals,a$CorlmCoru[t1[i]:t2[i],][,2],type='o',col='black',ylim=c(-0.25,0.25),xlab='Spearman',ylab=expression(cor[l] - cor[u]),cex.lab=2,cex.axis=1.5)
      arrows(spearvals, a$CorlmCoru[t1[i]:t2[i],][,1],spearvals,a$CorlmCoru[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='black')
      lines(range(spearvals),c(0,0),type='l',lty='dotted',col='black')
      par(new=T)
      plot(spearvals,t_CorlmCoru,pch=2,axes=F, xlab=NA, ylab=NA, col='black',ylim=c(0,1))
      lines(range(spearvals),c(0.05,0.05),type='l',lty='dashed',col='black')
      axis(side=4,col='black',col.axis="black",cex.axis=1.5)
      mtext(side = 4, line = 1.7, 'p values', col='black', cex=1.3)
      
      plot(spearvals,a$PlmPu[t1[i]:t2[i],][,2],type='o',col='black',ylim=c(-0.25,0.25),xlab='Spearman',ylab=expression(P[l] - P[u]),cex.lab=2,cex.axis=1.5) 
      arrows(spearvals,a$PlmPu[t1[i]:t2[i],][,1],spearvals,a$PlmPu[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='black')
      lines(range(spearvals),c(0,0),type='l',lty='dotted',col='black')  
      par(new=T)
      plot(spearvals,t_PlmPu,pch=2,axes=F, xlab=NA, ylab=NA, col='black',ylim=c(0,1))
      lines(range(spearvals),c(0.05,0.05),type='l',lty='dashed',col='black')
      axis(side=4,col='black',col.axis="black",cex.axis=1.5)
      mtext(side = 4, line = 1.7, 'p values', col='black', cex=1.3)
      
      
      plot(spearvals,a$D2umD2l[t1[i]:t2[i],][,2],type='o',col='black',ylim=c(-0.05,0.05),xlab='Spearman',ylab=expression(D[u]^2 - D[l]^2),cex.lab=2,cex.axis=1.5) 
      arrows(spearvals,a$D2umD2l[t1[i]:t2[i],][,1],spearvals,a$D2umD2l[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='black')
      lines(range(spearvals),c(0,0),type='l',lty='dotted',col='black')
      par(new=T)
      plot(spearvals,t_D2umD2l,pch=2,axes=F, xlab=NA, ylab=NA, col='black',ylim=c(0,1))
      lines(range(spearvals),c(0.05,0.05),type='l',lty='dashed',col='black')
      axis(side=4,col='black',col.axis="black",cex.axis=1.5)
      mtext(side = 4, line = 1.7, 'p values', col='black', cex=1.3)
      
    #}
    
  }
  
  op2<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("top", "1st column : Clayton, 2nd column : Normal, 3rd column : Frank, 4th column : Survival Clayton", 
         cex = 2, xpd = TRUE, horiz = TRUE, inset = c(0,0), 
         bty = "n") 
  par(op2)
  dev.off()
}

#-------------------------------------------------------------------------------------------------------
#                                              CODE ENDS HERE
#-------------------------------------------------------------------------------------------------------




