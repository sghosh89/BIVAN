#---------------------------------------------------------------
source("CopulaFunctions.R")
source("vivj_matrix.R")
source("good_loclist.R")
library(ncf)
#---------------------------------------------------------------
#Processing function for a copula approach to synchrony
#
#Args
#m : output from vivj_matrix.R
#
#Output - A list with these elements
#ranks        A dataframe with 2 columns, one for 
#               ds1 and one for ds2, corresponding to samples from the copula.
#spear        Spearman correlation (single number)
#kend         Kendall correlation (single number)
#Corl,Coru        covariance based statistics
#Pl,Pu      statistics show how the points within a copula scattered from right diagonal of the box within a pentagon
#Tl,Tu      statistics show how the points within a copula scattered along right diagonal of the box within a triangle
#D2l,D2u      measures the squared distance between the points from the right diagonal within a copula 

copsync<-function(m){

  vi<-m[,1]
  vj<-m[,2]
  
  #get mean and variance
  vi_mean<-mean(vi)
  vj_mean<-mean(vj)
  var_vi<-var(vi)
  var_vj<-var(vj)
  
  if (length(vi)>0){
    #get spear
    spear<-cor(vi,vj, method ="pearson") 
    
    #get kend
    kend<-cor(vi, vj, method ="kendall") 
    
    #----------------------------------------------------STATISTICS :2 ---------------------------------------------------------
    #get Cl, Cu   (covariance based new stat)
    stat2<-CorlCoru(vi,vj)
    Corl<-stat2[1]
    Coru<-stat2[2]
    
    #---------------------------------------------------------- STATISTICS : 4 -----------------------------------------------------------------
    #   get New statistics : Pl and Pu   # distance from right diagonal in lower and upper triangle based stat
    stat4<-PlPu(vi,vj)
    Sl_Su_Si_P<-stat4[[1]]
    Pl<-stat4[[2]]
    Pu<-stat4[[3]]
    
    #--------------------------------------------------------- STATISTICS : 6 -----------------------------------------------------------------------
    # get Rl : average of squared distance of points from the right diagonal of the box for lower triangle
    # get Ru : average of squared distance of points from the right diagonal of the box for upper triangle
    stat6<-D2lD2u(vi,vj)
    D2l<-stat6[1]
    D2u<-stat6[2]
    
  }else{
    
    spear<-NA
    kend<-NA
    Corl<-NA
    Coru<-NA
    Sl_Su_Si_P<-NA
    Pl<-NA
    Pu<-NA
    D2l<-NA
    D2u<-NA
    
  }
  
  return(list(ranks=data.frame(Rki=vi,Rkj=vj),
              spear=spear,kend=kend,
              Corl=Corl,Coru=Coru,
              Sl_Su_Si_P=Sl_Su_Si_P,Pl=Pl,Pu=Pu,
              D2l=D2l,D2u=D2u))
}

#---------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------
#This function takes a matrix, computes the mean of the non-diagonal 
#entries, and then resamples rows and columns to get confidence
#intervals for that mean. Meant to be applied to matrices of pairwise
#comparisons between sampling locations.
#
#1) square matrix M (the matrix, possible with non-NA diagonal entries)
#2) numresamp (the number off resamplings to do)
#3) prob (a vector of quantiles you want, default c(0.025,0.975))
#4) ploton  - produce the histogram plot or not? Default "FALSE"
#
#Returns: a list with these names elements
#1) datmean - The mean of the matrix except for the diagonals
#2) quantiles - as specified by prob
#3) fracgt0 - fraction of resampling-based values greater than 0
#4) fraclt0 - (will just be 1-fracgt0)
#
#Also generates the histogram plot if ploton=T

resampmn<-function(M,numresamp=10000,prob=c(.025,0.975),ploton=F)
{
  #get datmean
  diag(M)<-NA
  datmean<-mean(M,na.rm=T)
  
  #do the resampling
  resamp_stat<-c()
  for(j in 1:numresamp){
    resamp_rows<-sample(1:nrow(M),replace=T)
    M1<-M[resamp_rows,resamp_rows]
    r1<-mean(M1,na.rm=T)
    resamp_stat<-c(resamp_stat,r1)
    quantiles<-quantile(resamp_stat,prob=c(0.025,0.975))
    q2.5<-quantiles[1]
    q97.5<-quantiles[2]
  }
  
  #get fracgt0 : fraction of resampling-based values greater than 0
  fracgt0<-length(which(resamp_stat>0))/length(resamp_stat)
  
  #get fraclt0 : fraction of resampling-based values greater than 0
  #  fraclt0<-length(which(resamp_stat<0))/length(resamp_stat)
  
  
  #generate the plot
  if (ploton){
    hist(resamp_stat,breaks=100)
    abline(v=datmean,col="red",lwd=2)
    abline(v=q2.5,col="green4",lwd=1)
    abline(v=q97.5,col="blue",lwd=1)  
    mtext(paste0("datmean=",round(datmean,4),", q2.5=",round(q2.5,4),", q97.5=",round(q97.5,4))) 
  }
  
  #  return(list(datmean=datmean,quantiles=quantiles,fracgt0=fracgt0,fraclt0=fraclt0))
  return(list(datmean=datmean,quantiles=quantiles,fracgt0=fracgt0))
}

#---------------------------------------------------------------------------------------------------------------------------
#Calling the above on all pairs of several time series and plotting
#results and returning all stats. All this function does is call the 
#above function on all pairs of time series, organize the results, 
#and make plots.
#
#Args
#d_allsp            A list of data frames, each with columns Year and Dat
#               The years are assumed to be sequential and all included,
#               though there may be NAs in Dat and the years may not
#               be all the same for ds1 and ds2.
#sp             species number
#lats, longs  Vectors of latitudes and longitudes for the locations from
#               which data in d were gathered, both length equal to length(d)
#pfname       Filename (without extension) prepended to plot files saved.
#
#Output - A list with these elements
#D        A matrix of geographic distances between sampling locations
#spear    A matrix of spearman results, length(d) by length(d)
#kend     A matrix of kendall results, length(d) by length(d)
#Corl       A matrix of Cl results, length(d) by length(d)
#Coru       A matrix of Cu results, length(d) by length(d)
#Pl    A matrix of Shy_lt results, length(d) by length(d)
#Pu   A matrix of Shy_ut results, length(d) by length(d)

#D2l    A matrix of R_l results, length(d) by length(d)
#D2u    A matrix of R_u results, length(d) by length(d)
#numericdf    A dataframe containing all the statistical data

multcall<-function(d_allsp,sp,lats,longs,pfname,good_loc){
  
  d<-d_allsp[[sp]]
  lenloc<-length(good_loc)
  
  #D<-matrix(NA,lenloc,lenloc)
  D<-gcdist(longs[good_loc],lats[good_loc]) # This good_loc by good_loc matrix with distance btw any two selected among all locations
  colnames(D) <- paste("loc",good_loc, sep="")
  rownames(D) <- paste("loc",good_loc, sep="")
  
  #first initialize result receptacles for the output
  spear<-matrix(NA,lenloc,lenloc)
  colnames(spear) <- colnames(D)
  rownames(spear) <-rownames(D)
  
  kend<-matrix(NA,lenloc,lenloc)
  colnames(kend) <- colnames(D)
  rownames(kend) <-rownames(D)
  
  Corl<-matrix(NA,lenloc,lenloc)
  colnames(Corl) <- colnames(D)
  rownames(Corl) <-rownames(D)
  
  Coru<-matrix(NA,lenloc,lenloc)
  colnames(Coru) <- colnames(D)
  rownames(Coru) <-rownames(D)
  
  Pl<-matrix(NA,lenloc,lenloc)
  colnames(Pl) <- colnames(D)
  rownames(Pl) <-rownames(D)
  
  Pu<-matrix(NA,lenloc,lenloc)
  colnames(Pu) <- colnames(D)
  rownames(Pu) <-rownames(D)
  
  D2l<-matrix(NA,lenloc,lenloc)
  colnames(D2l) <- colnames(D)
  rownames(D2l) <-rownames(D)
  
  D2u<-matrix(NA,lenloc,lenloc)
  colnames(D2u) <- colnames(D)
  rownames(D2u) <-rownames(D)
  
  
  #------------------- PLOT :  copula_for_all_location pair ----------------
  pdf(paste(pfname,"_AllCops.pdf",sep=""),width=6*lenloc, height=6*lenloc)
  op<-par(mfrow=c(lenloc,lenloc),mar=c(3,3,3,3), mgp=c(1.5,0.5,0))
  
  for (ii in c(1:lenloc)){
    for (jj in c(1:lenloc)){
      #compute results
      i<-good_loc[ii]
      j<-good_loc[jj]
      #cat("i,j",i,j,"\n")
      
      #D[ii,jj]<-gcdist(longs[i],lats[i],longs[j],lats[j]) old way of gcdist calling from ncf-version 1.1 to 1.7
   
      m<-vivj_matrix(d_allsp,sp,i,j)
      thisres<-copsync(m)
      
      spear[ii,jj]<-thisres$spear
      kend[ii,jj]<-thisres$kend
      
      Corl[ii,jj]<-thisres$Corl
      Coru[ii,jj]<-thisres$Coru
      
      Pl[ii,jj]<-thisres$Pl
      Pu[ii,jj]<-thisres$Pu
      
      D2l[ii,jj]<-thisres$D2l
      D2u[ii,jj]<-thisres$D2u
      
      
      plot(thisres$ranks$Rki,thisres$ranks$Rkj,type='p',col=rgb(0,0,0,.2),pch=19,xlim=c(0,1),ylim=c(0,1),xlab=expression(u[i]),ylab=expression(v[j]),cex.lab=2)
      mtext(paste0("[ i, j ] ="," [",i,",",j,"] ", ","," n=",dim(thisres$ranks)[1],", D=",round(D[ii,jj],2)),side = 3, line=0.15, adj=0.5, col="red")
    }
  }
  
  par(op)
  dev.off()
  
  
  #------------------- PLOT :  (Sl, Su and Si)_P for_all_location pair ----------------
  pdf(paste(pfname,"_Sl_Su_Si_P.pdf",sep=""),width=6*lenloc, height=6*lenloc)
  op<-par(mfrow=c(lenloc,lenloc),mar=c(3,3,3,3), mgp=c(1.5,0.5,0))
  
  for (ii in 1:lenloc){
    for (jj in 1:lenloc){
      
      i<-good_loc[ii]
      j<-good_loc[jj]
      m<-vivj_matrix(d_allsp,sp,i,j)
      thisres<-copsync(m)
      
      if(is.na(thisres$Sl_Su_Si_P$Sl_P[1]) == F){    
        plot(thisres$Sl_Su_Si_P$dist_Sl_P,thisres$Sl_Su_Si_P$Sl_P,type='l',col="red",ylim=c(0,1),xlab=" ",ylab=" ")
        lines(thisres$Sl_Su_Si_P$dist_Su_P,thisres$Sl_Su_Si_P$Su_P,type='l',col="blue",ylim=c(0,1))
        lines(thisres$Sl_Su_Si_P$dist_Si_P,thisres$Sl_Su_Si_P$Si_P,type='l',lty="dashed",col="green4",ylim=c(0,1))
        mtext(paste0("[ i, j ] ="," [",i,",",j,"], "," red : Sl_P, "," blue : Su_P, ", " green : Si_P"),side = 3, line=0.15, adj=0.5, col="black")
      }else{
        plot(0,0,type='p',col="white",cex=0,xlim=c(0,sqrt(2)/2),ylim=c(0,1),xlab=" ",ylab=" ")
        mtext(paste0("[ i, j ] ="," [",i,",",j,"]," ),side = 3, line=0.15, adj=0.5, col="black")
        mtext(paste0("NA" ),side = 3, line=-22.5, adj=0.5, col="black", cex=4)
      }
    }
  }
  
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------
  numericdf<-data.frame(Stat=c('spear','kend','Corl','Coru','Pl','Pu','D2l','D2u','Corl-Coru','Pl-Pu','D2u-D2l'),
                        mnvalue=NA,Lower95CI=NA,Upper95CI=NA,fracgt0=NA)
  
  #---------------------- PLOT :  spear_vs_D -------------------------------
  pdf(paste(pfname,"_Spearman_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,spear,xlab="D",ylab="spear",col=rgb(0.5,0,0,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(spear)
  #***DAN: fill in one line of numericdf here, the line for spear
  i_name<-which(numericdf$Stat=="spear")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<spear>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  kend_vs_D ---------------------------------
  pdf(paste(pfname,"_Kendall_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,kend,xlab="D",ylab="kend",col=rgb(0,0,1,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(kend)
  i_name<-which(numericdf$Stat=="kend")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<kend>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  
  #--------------------- PLOT :  Corl_vs_D ---------------------------------
  pdf(paste(pfname,"_Corl_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,Corl,xlab="D",ylab=expression(cor[l]),col=rgb(0.5,0,0,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(Corl)
  i_name<-which(numericdf$Stat=="Corl")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<Corl>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  Coru_vs_D ---------------------------------
  pdf(paste(pfname,"_Coru_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,Coru,xlab="D",ylab=expression(cor[u]),col=rgb(0,0.5,0,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(Coru)
  i_name<-which(numericdf$Stat=="Coru")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<Coru>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  Corl-Coru_vs_D ---------------------------------
  pdf(paste(pfname,"_Corl-Coru_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,Corl-Coru,xlab="D",ylab=expression(cor[l]-cor[u]),col=rgb(0,0,0.5,.2),pch=19,ylim=c(-1,1),cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(Corl-Coru)
  i_name<-which(numericdf$Stat=="Corl-Coru")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<Corl-Coru>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------

  #--------------------- PLOT :  Pl_vs_D ---------------------------------
  pdf(paste(pfname,"_Pl_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,Pl,xlab="D",ylab=expression(P[l]),col=rgb(0.5,0,0,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(Pl)
  i_name<-which(numericdf$Stat=="Pl")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<Pl>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  Pu_vs_D ---------------------------------
  pdf(paste(pfname,"_Pu_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,Pu,xlab="D",ylab=expression(P[u]),col=rgb(0,0,0.5,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(Pu)
  i_name<-which(numericdf$Stat=="Pu")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<Pu>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  Pl-Pu_vs_D ---------------------------------
  pdf(paste(pfname,"_Pl-Pu_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,Pl-Pu,xlab="D",ylab=expression(P[l]-P[u]),col=rgb(0,0.5,0,.2),pch=19,ylim=c(-1,1),cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(Pl-Pu)
  i_name<-which(numericdf$Stat=="Pl-Pu")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<Pl-Pu>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  D2l_vs_D ---------------------------------
  pdf(paste(pfname,"_D2l_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,D2l,xlab="D",ylab=expression(D[l]^2),col=rgb(0.5,0,0,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(D2l)
  i_name<-which(numericdf$Stat=="D2l")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<D2l>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  D2u_vs_D ---------------------------------
  pdf(paste(pfname,"_D2u_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,D2u,xlab="D",ylab=expression(D[u]^2),col=rgb(0,0.5,0,.2),pch=19,cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(D2u)
  i_name<-which(numericdf$Stat=="D2u")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<D2u>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------
  
  #--------------------- PLOT :  D2u-D2l_vs_D ---------------------------------
  pdf(paste(pfname,"_D2u-D2l_vs_D.pdf",sep=""),width=8, height=8)
  op<-par(mgp=c(3.5,1,0),mar=c(7,7,1,1))
  plot(D,D2u-D2l,xlab="D",ylab=expression(D[u]^2 - D[l]^2),col=rgb(0,0,0.5,.2),pch=19,ylim=c(-0.3,0.3),cex.lab=3,cex.axis=2)
  lines(range(D),c(0,0),type='l',lty='dashed')
  result<-resampmn(D2u-D2l)
  i_name<-which(numericdf$Stat=="D2u-D2l")
  numericdf$mnvalue[i_name]<-result$datmean
  numericdf$Lower95CI[i_name]<-result$quantiles[1]
  numericdf$Upper95CI[i_name]<-result$quantiles[2]
  numericdf$fracgt0[i_name]<-result$fracgt0
  #mtext(paste0("<D2u-D2l>=",round(result$datmean,4),", q2.5=",round(result$quantiles[1],4),", q97.5=",round(result$quantiles[2],4),", fracgt0=",result$fracgt0))
  par(op)
  dev.off()
  #-------------------------------------------------------------------------

  return(list(D=D,
              spear=spear,kend=kend,
              Corl=Corl,Coru=Coru,
              Pl=Pl,Pu=Pu,
              D2l=D2l,D2u=D2u,
              numericdf=numericdf))
}

#----------------------------------------------------------------------------------------------------------------------------
#                                                       CODE ENDS HERE
#-----------------------------------------------------------------------------------------------------------------------------









