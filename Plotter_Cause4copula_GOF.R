#------------------------------------------------------------------------------------------------
# This function examines the model selection based goodness of fit results between noise copula input
# and pop. copula output (mechanism for non-N copula)
#-----------------------------------------------------------------------------------------------
# This function will generate a plot for C, SC, G and SG copula (or any other as the case may be) : 
#  parameter of noise copula vs. correln coef (spearman/kendall)
# BS : number of bootstrapps used for BiCopGOFTest

Plotter_Cause4copula_GOF<-function(N,fcode,method,nsd,num_keep_last,BS,params,p0,model){
  
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
    s<-GetNoise(N=N,fcode=fcode,corcoef=corcoef,nsd=nsd,method=method,ploton=F)
    # generate pop_copula
    s2<-Simulator_Cause4copula(params=params,p0=p0,noise=s$noise_q,model=model)
    
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
