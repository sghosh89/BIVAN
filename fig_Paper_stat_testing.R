#--------------------------------------------
source("./TestStats_multivar_dataset.R")
#--------------------------------------------

# This function make stat_testing plot for a specified copula with a specific statistic at a time
# It is actually the simpler version of "plotter_stat_testing" function from TestStats_multivar_dataset.R

call_fig_Paper_stat_testing<-function(resultsloc,statname,copname,reslist,numsims=100){
  
  pdf(paste(resultsloc,statname,"_",copname,".pdf",sep=""),width=4, height=3)
  a<-M1mM2(reslist=res_pt35_rankT,numsims)
  spearvals<-seq(from=0,to=0.9,by=0.1)
  callfn<-coplist_with_params(spearvals=spearvals)
  b<-callfn$paramlist
  
  t1<-c(1,length(b[[1]])+1,length(b[[1]])+length(b[[2]])+1,length(b[[1]])+length(b[[2]])+length(b[[3]])+1)
  t2<-c(length(b[[1]]),length(b[[1]])+length(b[[2]]),length(b[[1]])+length(b[[2]])+length(b[[3]]),length(b[[1]])+length(b[[2]])+length(b[[3]])+length(b[[4]]))
  
  
  if(copname=="C"){
    i<-1
  }else if(copname=="N"){
    i<-2
  }else if(copname=="F"){
    i<-3
  }else if(copname=="SC"){
    i<-4
  }else{
    i<-NA
    warning("specify copname",immediate.=T,call.=T)
  }
  
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
  
  op<-par(mar=c(3.5,4,1,3.5),mgp=c(2.2,1,0))
  
  if(statname=="CorlmCoru"){
    
    plot(spearvals,a$CorlmCoru[t1[i]:t2[i],][,2],type='b',col='orange',
         ylim=c(-0.25,0.25),xlab='Spearman',ylab=expression(cor[l] - cor[u]),cex.lab=1.2)
    arrows(spearvals, a$CorlmCoru[t1[i]:t2[i],][,1],spearvals,a$CorlmCoru[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='green')
    lines(range(spearvals),c(0,0),type='l',lty='dashed',col='deepskyblue1')
    par(new=T)
    plot(spearvals,t_CorlmCoru,pch=16,axes=F, xlab=NA, ylab=NA, col='purple2',ylim=c(0,1))
    lines(range(spearvals),c(0.05,0.05),type='l',lty='dashed',col='purple2')
    axis(side=4,col='purple2',col.axis="purple2")
    mtext(side = 4, line = 1.7, 'p values', col='purple2')
    
  }else if(statname=="PlmPu"){
    
    plot(spearvals,a$PlmPu[t1[i]:t2[i],][,2],type='b',col='black',
         ylim=c(-0.25,0.25),xlab='Spearman',ylab=expression(P[l] - P[u]),cex.lab=1.2) 
    arrows(spearvals,a$PlmPu[t1[i]:t2[i],][,1],spearvals,a$PlmPu[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='green')
    lines(range(spearvals),c(0,0),type='l',lty='dashed',col='deepskyblue1')  
    par(new=T)
    plot(spearvals,t_PlmPu,pch=16,axes=F, xlab=NA, ylab=NA, col='purple2',ylim=c(0,1))
    lines(range(spearvals),c(0.05,0.05),type='l',lty='dashed',col='purple2')
    axis(side=4,col='purple2',col.axis="purple2")
    mtext(side = 4, line = 1.7, 'p values', col='purple2')
    
  }else if(statname=="D2umD2l"){
    
    plot(spearvals,a$D2umD2l[t1[i]:t2[i],][,2],type='b',col='red',
         ylim=c(-0.05,0.05),xlab='Spearman',ylab=expression(D[u]^2 - D[l]^2),cex.lab=1.2) 
    arrows(spearvals,a$D2umD2l[t1[i]:t2[i],][,1],spearvals,a$D2umD2l[t1[i]:t2[i],][,3], length=0.02, angle=90, code=3, col='green')
    lines(range(spearvals),c(0,0),type='l',lty='dashed',col='deepskyblue1')
    par(new=T)
    plot(spearvals,t_D2umD2l,pch=16,axes=F, xlab=NA, ylab=NA, col='purple2',ylim=c(0,1))
    lines(range(spearvals),c(0.05,0.05),type='l',lty='dashed',col='purple2')
    axis(side=4,col='purple2',col.axis="purple2")
    mtext(side = 4, line = 1.7, 'p values', col='purple2')
    
  }else{
    plot(NA,NA,xlim=c(0,1),ylim=c(0,1))
    warning("specify statname",immediate.=T,call.=T)
  }
  par(op)
  
  dev.off()
  
}
#-------------------------------------------------------------------------------

