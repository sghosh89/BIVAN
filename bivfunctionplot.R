source("CopulaFunctions_flexible.R")
source("copsurrog2d.R")
#-----------------------------------------------------------
#utility function for calculating all the stats we want
calcstats<-function(v,f,nm,numbin){
  
  #calculate the stats
  res<-c()
  bingap<-1/numbin
  x<-seq(from=0,to=1,by=bingap)
  for (i in 1:(length(x)-1)){
    lbd<-x[i]
    #print(lbd)
    res<-c(res,f(v[,1],v[,2],lbd,lbd+bingap))
  }
  
  x1<-head(x,n=-1)
  x2<-tail(x,n=-1)
  #name the result vector
  h<-paste(nm,x1,"to",
           x2,sep='')
  h<-gsub(".","p",h,fixed=T)
  names(res)<-h
  
  return(res)
}
#-----------------------------------------------------
makeSurrog<-function(v,numsurrog=1000){
  cop<-normalCopula(.5,2)
  #numsurrog<-1000
  surv_K<-copsurrog2d(v,cop,"kendall",numsurrog) #Do the kendall surrogates
  surv_S<-copsurrog2d(v,cop,"spearman",numsurrog) #Do the spearman surrogates
  return(list(surv_K=surv_K,
              surv_S=surv_S))
}
#makeSurrog(v=v_BMR)
#------------------------------------------------------------------
#utility function for calculating fractions of surrogates 
#with a stat smaller than that of data, i.e., fraction of 
#surrogates for which the stat value on the data is bigger
fracwork<-function(dvals,surrvals){
  frac<-NA*numeric(length(dvals))
  names(frac)<-names(dvals)
  for (counter in 1:length(frac)){
    frac[counter]<-sum(surrvals[counter,]<dvals[counter])
  }
  return(frac)
}
#corstats_frac_K<-fracwork(corstats_d,corstats_K)

#--------------------------------------------------------------------
#now work with the P stats
Pbds_wrap<-function(vi,vj,lb,ub){
  return(Pbds(vi,vj,lb,ub)$abs_res)
}
#----------------------------------------------------
fracplot<-function(x,fracs,ylims,numsurrog=1000){
  for (counter in 1:length(fracs)){
    ptxt<-''
    if ((fracs[counter]>.975*numsurrog)==T){
      ptxt<-paste(">",fracs[counter],sep='')
    } 
    if ((fracs[counter]<.025*numsurrog)==T){
      ptxt<-paste("<",numsurrog-fracs[counter],sep='')
    }
    yht<-ylims[2]-.1*diff(ylims)
    text(x[counter],yht,ptxt,adj=c(0.5,.5),cex=0.5,srt=90)
  }
}
#-----------------------------
#This function will do the same thing as of fracplot
myfracplot<-function(x,stats_d,surrog_statsq,stats_frac,ylims,numsurrog=1000){
  for (counter in 1:length(x)){
    ptxt<-''
    if(stats_d[counter]>surrog_statsq[3,][counter]){
      ptxt<-paste(">",stats_frac[counter],sep='')
    }
    if(stats_d[counter]<surrog_statsq[2,][counter]){
      ptxt<-paste("<",numsurrog-stats_frac[counter],sep='')
    }
    #ptxt2<-stats_frac[counter]
    
    yht<-ylims[2]-.1*diff(ylims)
    text(x[counter],yht,ptxt,adj=c(0.5,.5),cex=0.5,srt=90)
    #text(x[counter],yht,ptxt2,adj=c(0.5,.5),cex=0.5,srt=90)
  }
}
#----------------------------------------------------------------
vlineplot<-function(x,ylims){
  for (counter in 1:length(x)){
    lines(rep(x[counter]-(x[2]-x[1])/2,2),ylims,type='l',lty='dotted')
  }
  lines(rep(x[length(x)]+(x[2]-x[1])/2,2),ylims,type='l',lty='dotted')
}
#---------------------------------------------------------------------
# This is the plotter function as well as it returns the stats as a list of 6
bivfunctionplot<-function(v,resloc,nametag,numbin,numsurrog=1000){
  
  temp<-makeSurrog(v=v)
  surv_K<-temp$surv_K
  surv_S<-temp$surv_S
  
  #Cor stats
  corstats_d<-calcstats(v=v,f=Corbds,nm="Cor",numbin=numbin)
  corstats_K<-apply(FUN=calcstats,X=surv_K,MARGIN=3,f=Corbds,nm="Cor",numbin=numbin)
  corstats_frac_K<-fracwork(corstats_d,corstats_K)
  
  corlmcoru_d<-unname(corstats_d[1]-corstats_d[length(corstats_d)])
  corlmcoru_K<-corstats_K[1,]-corstats_K[length(corstats_d),]
  corlmcoru_frac_K<-sum(corlmcoru_K<corlmcoru_d)
  corstats_Kq<-apply(FUN=quantile,X=corstats_K,MARGIN=1,prob=c(.005,0.025,.975,.995))
  
  corstats_S<-apply(FUN=calcstats,X=surv_S,MARGIN=3,f=Corbds,nm="Cor",numbin=numbin)
  corstats_frac_S<-fracwork(corstats_d,corstats_S)
  
  corlmcoru_S<-corstats_S[1,]-corstats_S[length(corstats_d),]
  corlmcoru_frac_S<-sum(corlmcoru_S<corlmcoru_d)
  corstats_Sq<-apply(FUN=quantile,X=corstats_S,MARGIN=1,prob=c(.005,0.025,.975,.995))

  #P stats
  Pstats_d<-calcstats(v=v,f=Pbds_wrap,nm="P",numbin=numbin)
  Pstats_K<-apply(FUN=calcstats,X=surv_K,MARGIN=3,f=Pbds_wrap,nm="P",numbin=numbin)
  Pstats_frac_K<-fracwork(Pstats_d,Pstats_K)
  
  PlmPu_d<-unname(Pstats_d[1]-Pstats_d[length(Pstats_d)])
  PlmPu_K<-Pstats_K[1,]-Pstats_K[length(Pstats_d),]
  PlmPu_frac_K<-sum(PlmPu_K<PlmPu_d)
  Pstats_Kq<-apply(FUN=quantile,X=Pstats_K,MARGIN=1,prob=c(.005,0.025,.975,.995))
  
  Pstats_S<-apply(FUN=calcstats,X=surv_S,MARGIN=3,f=Pbds_wrap,nm="P",numbin=numbin)
  Pstats_frac_S<-fracwork(Pstats_d,Pstats_S)
  
  PlmPu_S<-Pstats_S[1,]-Pstats_S[length(Pstats_d),]
  PlmPu_frac_S<-sum(PlmPu_S<PlmPu_d)
  Pstats_Sq<-apply(FUN=quantile,X=Pstats_S,MARGIN=1,prob=c(.005,0.025,.975,.995))
  
  #-------------------------------Ranking--------------------------
  Rank_Pl_K<-unname((numsurrog-Pstats_frac_K)[1])
  Rank_Pl_S<-unname((numsurrog-Pstats_frac_S)[1])
  Rank_Pu_K<-unname((numsurrog-Pstats_frac_K)[length(Pstats_frac_K)])
  Rank_Pu_S<-unname((numsurrog-Pstats_frac_S)[length(Pstats_frac_K)])
  Rank_PlmPu_S<-numsurrog-PlmPu_frac_S
  Rank_PlmPu_K<-numsurrog-PlmPu_frac_K
  #---------------------------------------------------------------------
  
  #now work with the D2 stats
  D2stats_d<-calcstats(v=v,f=D2bds,nm="Dtwo",numbin=numbin)
  D2stats_K<-apply(FUN=calcstats,X=surv_K,MARGIN=3,f=D2bds,nm="Dtwo",numbin=numbin)
  D2stats_frac_K<-fracwork(D2stats_d,D2stats_K)
  
  D2umD2l_d<-D2stats_d[length(D2stats_d)]-D2stats_d[1]
  D2umD2l_K<-D2stats_K[length(D2stats_d),]-D2stats_K[1,]
  D2umD2l_frac_K<-sum(D2umD2l_K<D2umD2l_d)
  D2stats_Kq<-apply(FUN=quantile,X=D2stats_K,MARGIN=1,prob=c(.005,0.025,.975,.995))
  
  D2stats_S<-apply(FUN=calcstats,X=surv_S,MARGIN=3,f=D2bds,nm="Dtwo",numbin=numbin)
  D2stats_frac_S<-fracwork(D2stats_d,D2stats_S)
  
  D2umD2l_S<-D2stats_S[length(D2stats_d),]-D2stats_S[1,]
  D2umD2l_frac_S<-sum(D2umD2l_S<D2umD2l_d)
  D2stats_Sq<-apply(FUN=quantile,X=D2stats_S,MARGIN=1,prob=c(.005,0.025,.975,.995))
  
  # -----------------Plot---------------------------------------------------------------------
  bingap<-1/numbin
  x<-seq(from=0,to=1,by=bingap)
  x<-head(x,-1)+diff(x)/2
  #x<-seq(from=0.05,to=0.95,by=0.1)
  xlimits<-c(0,1)
  
  #plotting layout, units inches
  xht<-0.5
  ywd<-0.5
  titleht<-.25
  panht<-1.5
  panwd<-3
  gap<-.05
  totwd<-ywd+2*panwd+2*gap
  totht<-xht+3*panht+3*gap+titleht
  pdf(paste(resloc,nametag,"_bivfunctionplot.pdf",sep=""),width=totwd,height=totht)
  #pdf(file="./Results/stat_results/stat_soilCN/CNdata_NonparamTailDep.pdf",width=totwd,height=totht)
  
  #plot kendall results on the left panels, spearman on right,
  #cor in top panels, then P, then D2 on bottom panels
  
  #kendall, cor
  par(fig=c(ywd/totwd,
            (ywd+panwd)/totwd,
            (xht+2*panht+2*gap)/totht,
            (xht+3*panht+2*gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
  ylimits_cor<-range(corstats_d,corstats_S,corstats_K)#??? should it be range(corstats_d,corstats_S,corstats_K)???
  ylimits_cor[2]<-ylimits_cor[2]+.3*diff(ylimits_cor)
  plot(x,corstats_d,type='p',pch=3,col="red",xlim=xlimits,ylim=ylimits_cor,
       xaxt='n',cex=1.5)
  mtext(side=3,line=0,text="Kendall-preserving surrogates")
  mtext(side=2,line=1,text="Partial correlation")
  axis(side=1,labels=F)
  #segments(x0=x,y0=corstats_K[2,],x1=x,y1=corstats_K[3,],col="blue")
  #lines(x,corstats_K[1,],type='l',lty='dotted')
  points(x,corstats_Kq[2,],pch=4,col="blue") # low CI 0.025
  points(x,corstats_Kq[3,],pch=4,col="green") # up CI 0.975
  #lines(x,corstats_K[4,],type='l',lty='dotted')
  text(xlimits[1],ylimits_cor[1],labels='A',cex=1.5,adj=c(.5,0))
  #myfracplot(x=x,stats_d = corstats_d,surrog_statsq = corstats_Kq,stats_frac = corstats_frac_K,
  #          ylims = ylimits_cor,numsurrog = 1000)
  fracplot(x=x,corstats_frac_K,ylimits_cor)
  vlineplot(x,ylimits_cor)
  
  #spearman, cor
  par(fig=c((ywd+panwd+gap)/totwd,
            (ywd+2*panwd+gap)/totwd,
            (xht+2*panht+2*gap)/totht,
            (xht+3*panht+2*gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
  plot(x,corstats_d,type='p',pch=3,col="red",xlim=xlimits,ylim=ylimits_cor,
       xaxt='n',yaxt='n',cex=1.5)
  mtext(side=3,line=0,text="Spearman-preserving surrogates")
  axis(side=1,labels=F)
  axis(side=2,labels=F)
  #lines(x,corstats_S[1,],type='l',lty='dotted')
  points(x,corstats_Sq[2,],pch=4,col="blue")
  points(x,corstats_Sq[3,],pch=4,col="green")
  #lines(x,corstats_S[4,],type='l',lty='dotted')
  text(xlimits[1],ylimits_cor[1],labels='B',cex=1.5,adj=c(.5,0))
  #myfracplot(x=x,stats_d = corstats_d,surrog_statsq = corstats_Sq,stats_frac = corstats_frac_S,
  #           ylims = ylimits_cor,numsurrog = 1000)
  fracplot(x=x,corstats_frac_S,ylimits_cor)
  vlineplot(x,ylimits_cor)
  
  #kendall, P
  par(fig=c(ywd/totwd,
            (ywd+panwd)/totwd,
            (xht+panht+gap)/totht,
            (xht+2*panht+gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
  ylimits_P<-range(Pstats_d,Pstats_S,Pstats_K)
  ylimits_P[2]<-ylimits_P[2]+.3*diff(ylimits_P)
  plot(x,Pstats_d,type='p',pch=3,col="red",xlim=xlimits,ylim=ylimits_P,xaxt='n',cex=1.5)
  axis(side=1,labels=F)
  mtext(side=2,line=1,text="P")
  #lines(x,Pstats_K[1,],type='l',lty='dotted')
  points(x,Pstats_Kq[2,],pch=4,col="blue")
  points(x,Pstats_Kq[3,],pch=4,col="green")
  #lines(x,Pstats_K[4,],type='l',lty='dotted')
  text(xlimits[1],ylimits_P[1],labels='C',cex=1.5,adj=c(.5,0))
  #myfracplot(x=x,stats_d = Pstats_d,surrog_statsq = Pstats_Kq,stats_frac = Pstats_frac_K,
  #           ylims = ylimits_P,numsurrog = 1000)
  fracplot(x=x,Pstats_frac_K,ylimits_P)
  vlineplot(x,ylimits_P)
  
  #spearman, p
  par(fig=c((ywd+panwd+gap)/totwd,
            (ywd+2*panwd+gap)/totwd,
            (xht+panht+gap)/totht,
            (xht+2*panht+gap)/totht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
  plot(x,Pstats_d,type='p',pch=3,col="red",xlim=xlimits,ylim=ylimits_P,
       xaxt='n',yaxt='n',cex=1.5)
  axis(side=1,labels=F)
  axis(side=2,labels=F)
  #lines(x,Pstats_S[1,],type='l',lty='dotted')
  points(x,Pstats_Sq[2,],pch=4,col="blue")
  points(x,Pstats_Sq[3,],pch=4,col="green")
  #lines(x,Pstats_S[4,],type='l',lty='dotted')
  text(xlimits[1],ylimits_P[1],labels='D',cex=1.5,adj=c(.5,0))
  #myfracplot(x=x,stats_d = Pstats_d,surrog_statsq = Pstats_Sq,stats_frac = Pstats_frac_S,
  #           ylims = ylimits_P,numsurrog = 1000)
  fracplot(x=x,Pstats_frac_S,ylimits_P)
  vlineplot(x,ylimits_P)
  
  #kendall, D2
  par(fig=c(ywd/totwd,
            (ywd+panwd)/totwd,
            (xht)/totht,
            (xht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
  ylimits_D2<-range(D2stats_d,D2stats_S,D2stats_K)
  ylimits_D2[2]<-ylimits_D2[2]+.3*diff(ylimits_D2)
  plot(x,D2stats_d,type='p',pch=3,col="red",xlim=xlimits,ylim=ylimits_D2,cex=1.5)
  mtext(side=1,line=1,text="Diagonal slice")
  mtext(side=2,line=1,text=expression(D^{2}))
  axis(side=1,labels=F)
  #lines(x,D2stats_K[1,],type='l',lty='dotted')
  points(x,D2stats_Kq[2,],pch=4,col="blue")
  points(x,D2stats_Kq[3,],pch=4,col="green")
  #lines(x,D2stats_K[4,],type='l',lty='dotted')
  text(xlimits[1],ylimits_D2[1],labels='E',cex=1.5,adj=c(.5,0))
  #myfracplot(x=x,stats_d = D2stats_d,surrog_statsq = D2stats_Kq,stats_frac = D2stats_frac_K,
  #           ylims = ylimits_D2,numsurrog = 1000)
  fracplot(x=x,D2stats_frac_K,ylimits_D2)
  vlineplot(x,ylimits_D2)
  
  #spearman, D2
  par(fig=c((ywd+panwd+gap)/totwd,
            (ywd+2*panwd+gap)/totwd,
            (xht)/totht,
            (xht+panht)/totht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
  plot(x,D2stats_d,type='p',pch=3,col="red",xlim=xlimits,ylim=ylimits_D2,
       yaxt='n',cex=1.5)
  mtext(side=1,line=1,text="Diagonal slice")
  axis(side=1,labels=F)
  axis(side=2,labels=F)
  #lines(x,D2stats_S[1,],type='l',lty='dotted')
  points(x,D2stats_Sq[2,],pch=4,col="blue")
  points(x,D2stats_Sq[3,],pch=4,col="green")
  #lines(x,D2stats_S[4,],type='l',lty='dotted')
  text(xlimits[1],ylimits_D2[1],labels='F',cex=1.5,adj=c(.5,0))
  #myfracplot(x=x,stats_d = D2stats_d,surrog_statsq = D2stats_Sq,stats_frac = D2stats_frac_S,
  #           ylims = ylimits_D2,numsurrog = 1000)
  fracplot(x=x,D2stats_frac_S,ylimits_D2)
  vlineplot(x,ylimits_D2)
  
  dev.off()
  
  return(list(corlmcoru_frac_K=corlmcoru_frac_K,
              corlmcoru_frac_S=corlmcoru_frac_S,
              PlmPu_frac_K=PlmPu_frac_K,
              PlmPu_frac_S=PlmPu_frac_S,
              D2umD2l_frac_K=D2umD2l_frac_K,
              D2umD2l_frac_S=D2umD2l_frac_S,
              Rank_Pl_K=Rank_Pl_K,
              Rank_Pl_S=Rank_Pl_S,
              Rank_Pu_K=Rank_Pu_K,
              Rank_Pu_S=Rank_Pu_S,
              Rank_PlmPu_K=Rank_PlmPu_K,
              Rank_PlmPu_S=Rank_PlmPu_S))
  
}
#-------------------------------------------------------------------------
#set.seed(seed=101) 
#source("getcopula.R")
#d<-readRDS("Data/RaCA_soilorganicC_soiltotalN_stocks100cm.RDS")
#d<-d[,c("SOCstock100","TSNstock100")]
#v_CN<-getcopula(d=d,rankon=T,ploton=T) 
#xxx<-bivfunctionplot(v=v_CN,resloc="./Results/stat_results/stat_soilCN/",nametag="trial",numbin=10)








