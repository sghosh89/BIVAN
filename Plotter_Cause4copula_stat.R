#This is a plotter function to plot the outputs from Sim_Cause4copula_stat function written in Cause4copula.R
# Input :
#       res : a big list which is the output from Sim_Cause4copula_stat function

Plotter_Cause4copula_stat<-function(res){
  
  resloc<-res$resloc
  method<-res$method
  fcode<-res$fcode
  corcoef_list<-res$corcoef_list
  S_noise_mat<-res$S_noise_mat
  K_noise_mat<-res$K_noise_mat
  P_noise_mat<-res$P_noise_mat
  S_pop_mat<-res$S_pop_mat
  K_pop_mat<-res$K_pop_mat
  P_pop_mat<-res$P_pop_mat
  Corl_noise_mat<-res$Corl_noise_mat
  Corl_pop_mat<-res$Corl_pop_mat
  Coru_noise_mat<-res$Coru_noise_mat
  Coru_pop_mat<-res$Coru_pop_mat
  CorlmCoru_noise_mat<-res$CorlmCoru_noise_mat
  CorlmCoru_pop_mat<-res$CorlmCoru_pop_mat
  Pl_noise_mat<-res$Pl_noise_mat
  Pl_pop_mat<-res$Pl_pop_mat
  Pu_noise_mat<-res$Pu_noise_mat
  Pu_pop_mat<-res$Pu_pop_mat
  PlmPu_noise_mat<-res$PlmPu_noise_mat
  PlmPu_pop_mat<-res$PlmPu_pop_mat
  D2u_noise_mat<-res$D2u_noise_mat
  D2u_pop_mat<-res$D2u_pop_mat
  D2l_noise_mat<-res$D2l_noise_mat
  D2l_pop_mat<-res$D2l_pop_mat
  D2umD2l_noise_mat<-res$D2umD2l_noise_mat
  D2umD2l_pop_mat<-res$D2umD2l_pop_mat
  pval_S<-res$pval_S
  pval_K<-res$pval_K
  pval_P<-res$pval_P
  pval_Corl<-res$pval_Corl
  pval_Coru<-res$pval_Coru
  pval_CorlmCoru<-res$pval_CorlmCoru
  pval_Pl<-res$pval_Pl
  pval_Pu<-res$pval_Pu
  pval_PlmPu<-res$pval_PlmPu
  pval_D2u<-res$pval_D2u
  pval_D2l<-res$pval_D2l
  pval_D2umD2l<-res$pval_D2umD2l
  
  
  if(method=="spearman"){
    xlabel<-"Spearman's Rho"
  }else if(method=="kendall"){
    xlabel<-"Kendall's Tau"
  }else{
    warning("specify method",immediate.=T,call.=T)
  }
  
  # Plotting Spearman correlation btw noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Spearman_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,S_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab="Spearman",xlim=c(0,1),ylim=c(0,1),cex.lab=2,cex.axis=1.5)
  segments(corcoef_list,S_noise_mat[,1],corcoef_list,S_noise_mat[,3],col='dimgrey')
  bar_len<-0.02
  segments(corcoef_list-bar_len,S_noise_mat[,1],corcoef_list+bar_len,S_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,S_noise_mat[,3],corcoef_list+bar_len,S_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,S_noise_mat[,1],corcoef_list,S_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,S_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,S_pop_mat[,1],corcoef_list,S_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,S_pop_mat[,1],corcoef_list+bar_len,S_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,S_pop_mat[,3],corcoef_list+bar_len,S_pop_mat[,3],col='black')
  #arrows(corcoef_list,S_pop_mat[,1],corcoef_list,S_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  points(corcoef_list,pval_S,col="black",pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting Kendall correlation btw noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Kendall_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,K_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab="Kendall",xlim=c(0,1),ylim=c(0,1),cex.lab=2,cex.axis=1.5)
  segments(corcoef_list,K_noise_mat[,1],corcoef_list,K_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,K_noise_mat[,1],corcoef_list+bar_len,K_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,K_noise_mat[,3],corcoef_list+bar_len,K_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,K_noise_mat[,1],corcoef_list,K_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,K_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,K_pop_mat[,1],corcoef_list,K_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,K_pop_mat[,1],corcoef_list+bar_len,K_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,K_pop_mat[,3],corcoef_list+bar_len,K_pop_mat[,3],col='black')
  #arrows(corcoef_list,K_pop_mat[,1],corcoef_list,K_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  points(corcoef_list,pval_K,col="black",pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting Pearson correlation btw noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pearson_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,P_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab="Pearson",xlim=c(0,1),ylim=c(0,1),cex.lab=2,cex.axis=1.5)
  segments(corcoef_list,P_noise_mat[,1],corcoef_list,P_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,P_noise_mat[,1],corcoef_list+bar_len,P_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,P_noise_mat[,3],corcoef_list+bar_len,P_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,P_noise_mat[,1],corcoef_list,P_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,P_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,P_pop_mat[,1],corcoef_list,P_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,P_pop_mat[,1],corcoef_list+bar_len,P_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,P_pop_mat[,3],corcoef_list+bar_len,P_pop_mat[,3],col='black')
  #arrows(corcoef_list,P_pop_mat[,1],corcoef_list,P_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  points(corcoef_list,pval_P,col="black",pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting Corl of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Corl_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Corl_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("Cor"["l"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.02+max(Corl_noise_mat[,2],Corl_pop_mat[,2])))
  segments(corcoef_list,Corl_noise_mat[,1],corcoef_list,Corl_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,Corl_noise_mat[,1],corcoef_list+bar_len,Corl_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,Corl_noise_mat[,3],corcoef_list+bar_len,Corl_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,Corl_noise_mat[,1],corcoef_list,Corl_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,Corl_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,Corl_pop_mat[,1],corcoef_list,Corl_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,Corl_pop_mat[,1],corcoef_list+bar_len,Corl_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,Corl_pop_mat[,3],corcoef_list+bar_len,Corl_pop_mat[,3],col='black')
  #arrows(corcoef_list,Corl_pop_mat[,1],corcoef_list,Corl_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_Corl,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting Coru of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Coru_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Coru_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("Cor"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.02+max(Coru_noise_mat[,2],Coru_pop_mat[,2])))
  segments(corcoef_list,Coru_noise_mat[,1],corcoef_list,Coru_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,Coru_noise_mat[,1],corcoef_list+bar_len,Coru_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,Coru_noise_mat[,3],corcoef_list+bar_len,Coru_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,Coru_noise_mat[,1],corcoef_list,Coru_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,Coru_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,Coru_pop_mat[,1],corcoef_list,Coru_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,Coru_pop_mat[,1],corcoef_list+bar_len,Coru_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,Coru_pop_mat[,3],corcoef_list+bar_len,Coru_pop_mat[,3],col='black')
  #arrows(corcoef_list,Coru_pop_mat[,1],corcoef_list,Coru_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_Coru,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting Corl - Coru of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Corl-Coru_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,CorlmCoru_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("Cor"["l"]-"Cor"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(-0.15,0.15))
  segments(corcoef_list,CorlmCoru_noise_mat[,1],corcoef_list,CorlmCoru_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,CorlmCoru_noise_mat[,1],corcoef_list+bar_len,CorlmCoru_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,CorlmCoru_noise_mat[,3],corcoef_list+bar_len,CorlmCoru_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,CorlmCoru_noise_mat[,1],corcoef_list,CorlmCoru_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,CorlmCoru_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,CorlmCoru_pop_mat[,1],corcoef_list,CorlmCoru_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,CorlmCoru_pop_mat[,1],corcoef_list+bar_len,CorlmCoru_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,CorlmCoru_pop_mat[,3],corcoef_list+bar_len,CorlmCoru_pop_mat[,3],col='black')
  #arrows(corcoef_list,CorlmCoru_pop_mat[,1],corcoef_list,CorlmCoru_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  lines(range(0,1),c(0,0),type='l',lty='dotted',col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_CorlmCoru,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting correlation between Corl - Coru of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_scatter_CorlmCoru_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(4.5,6,2,2), mgp=c(3,0.5,0))
  ylim1<-range(CorlmCoru_pop_mat[,1],CorlmCoru_pop_mat[,3])
  xlim1<-range(CorlmCoru_noise_mat[,1],CorlmCoru_noise_mat[,3])
  b_len<-diff(xlim1)/50
  plot(CorlmCoru_noise_mat[,2],CorlmCoru_pop_mat[,2],cex=2,col="black",
       xlab=expression("Cor"["l,input"]-"Cor"["u,input"]),
       ylab=expression("Cor"["l,output"]-"Cor"["u,output"]),
       cex.lab=2,cex.axis=1.5,
       xlim=xlim1,ylim=ylim1)
  segments(CorlmCoru_noise_mat[,2],CorlmCoru_pop_mat[,1],CorlmCoru_noise_mat[,2],CorlmCoru_pop_mat[,3],col='black')
  segments(CorlmCoru_noise_mat[,2]-b_len,CorlmCoru_pop_mat[,1],CorlmCoru_noise_mat[,2]+b_len,CorlmCoru_pop_mat[,1],col='black')
  segments(CorlmCoru_noise_mat[,2]-b_len,CorlmCoru_pop_mat[,3],CorlmCoru_noise_mat[,2]+b_len,CorlmCoru_pop_mat[,3],col='black')
  
  segments(CorlmCoru_noise_mat[,1],CorlmCoru_pop_mat[,2],CorlmCoru_noise_mat[,3],CorlmCoru_pop_mat[,2],col='dimgrey')
  segments(CorlmCoru_noise_mat[,1],CorlmCoru_pop_mat[,2]-b_len,CorlmCoru_noise_mat[,1],CorlmCoru_pop_mat[,2]+b_len,col='dimgrey')
  segments(CorlmCoru_noise_mat[,3],CorlmCoru_pop_mat[,2]-b_len,CorlmCoru_noise_mat[,3],CorlmCoru_pop_mat[,2]+b_len,col='dimgrey')
  
  dat<-data.frame(x1=CorlmCoru_noise_mat[,2],y1=CorlmCoru_pop_mat[,2])
  mylm<-lm(y1~x1,data=dat)
  abline(mylm,col="black")
  lines(x=xlim1,y=ylim1,col="black",lty="dashed")
  c<-cor.test(dat$x1,dat$y1,method = "pearson",alternative = "t")
  mtext(paste0("Pearson = ",round(unname(c$estimate),3),", p = ",round(c$p.value,3),sep=""),cex=1.5,line=0.1)
  par(op)
  dev.off()
  
  # Plotting Pl of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pl_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Pl_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("P"["l"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.01+max(Pl_noise_mat[,2],Pl_pop_mat[,2])))
  segments(corcoef_list,Pl_noise_mat[,1],corcoef_list,Pl_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,Pl_noise_mat[,1],corcoef_list+bar_len,Pl_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,Pl_noise_mat[,3],corcoef_list+bar_len,Pl_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,Pl_noise_mat[,1],corcoef_list,Pl_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,Pl_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,Pl_pop_mat[,1],corcoef_list,Pl_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,Pl_pop_mat[,1],corcoef_list+bar_len,Pl_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,Pl_pop_mat[,3],corcoef_list+bar_len,Pl_pop_mat[,3],col='black')
  #arrows(corcoef_list,Pl_pop_mat[,1],corcoef_list,Pl_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_Pl,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  #op<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 1, 0, 0), new = TRUE)
  #plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #legend("top", c("noise copula","population copula", "p-value(paired t-test)"), col = c("dimgrey", "black", "black"),
  #       cex = 0.8, pch = c(16, 16, 1), xpd = TRUE, horiz = TRUE, inset = c(0,0),
  #       bty = "n") 
  #par(op)
  dev.off()
  
  # Plotting Pu of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pu_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,Pu_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("P"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.01+max(Pu_noise_mat[,2],Pu_pop_mat[,2])))
  segments(corcoef_list,Pu_noise_mat[,1],corcoef_list,Pu_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,Pu_noise_mat[,1],corcoef_list+bar_len,Pu_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,Pu_noise_mat[,3],corcoef_list+bar_len,Pu_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,Pu_noise_mat[,1],corcoef_list,Pu_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,Pu_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,Pu_pop_mat[,1],corcoef_list,Pu_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,Pu_pop_mat[,1],corcoef_list+bar_len,Pu_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,Pu_pop_mat[,3],corcoef_list+bar_len,Pu_pop_mat[,3],col='black')
  #arrows(corcoef_list,Pu_pop_mat[,1],corcoef_list,Pu_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_Pu,col="black",type="p",axes = FALSE,
       bty = "n", xlab = "", ylab = "",ylim=c(0,1),xlim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting Pl - Pu of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_Pl-Pu_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,PlmPu_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("P"["l"]-"P"["u"]),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(-0.04,0.04))
  segments(corcoef_list,PlmPu_noise_mat[,1],corcoef_list,PlmPu_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,PlmPu_noise_mat[,1],corcoef_list+bar_len,PlmPu_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,PlmPu_noise_mat[,3],corcoef_list+bar_len,PlmPu_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,PlmPu_noise_mat[,1],corcoef_list,PlmPu_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,PlmPu_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,PlmPu_pop_mat[,1],corcoef_list,PlmPu_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,PlmPu_pop_mat[,1],corcoef_list+bar_len,PlmPu_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,PlmPu_pop_mat[,3],corcoef_list+bar_len,PlmPu_pop_mat[,3],col='black')
  #arrows(corcoef_list,PlmPu_pop_mat[,1],corcoef_list,PlmPu_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  lines(range(0,1),c(0,0),type='l',lty='dotted',col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_PlmPu,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting correlation between Pl - Pu of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_scatter_PlmPu_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(4.5,6,2,2), mgp=c(3,0.5,0))
  ylim1<-range(PlmPu_pop_mat[,1],PlmPu_pop_mat[,3])
  xlim1<-range(PlmPu_noise_mat[,1],PlmPu_noise_mat[,3])
  
  b_len<-diff(xlim1)/50
  
  plot(PlmPu_noise_mat[,2],PlmPu_pop_mat[,2],cex=2,col="black",
       xlab=expression("P"["l,input"]-"P"["u,input"]),
       ylab=expression("P"["l,output"]-"P"["u,output"]),
       cex.lab=2,cex.axis=1.5,
       xlim=xlim1,ylim=ylim1)
  
  segments(PlmPu_noise_mat[,2],PlmPu_pop_mat[,1],PlmPu_noise_mat[,2],PlmPu_pop_mat[,3],col='black')
  segments(PlmPu_noise_mat[,2]-b_len,PlmPu_pop_mat[,1],PlmPu_noise_mat[,2]+b_len,PlmPu_pop_mat[,1],col='black')
  segments(PlmPu_noise_mat[,2]-b_len,PlmPu_pop_mat[,3],PlmPu_noise_mat[,2]+b_len,PlmPu_pop_mat[,3],col='black')
  
  segments(PlmPu_noise_mat[,1],PlmPu_pop_mat[,2],PlmPu_noise_mat[,3],PlmPu_pop_mat[,2],col='dimgrey')
  segments(PlmPu_noise_mat[,1],PlmPu_pop_mat[,2]-b_len,PlmPu_noise_mat[,1],PlmPu_pop_mat[,2]+b_len,col='dimgrey')
  segments(PlmPu_noise_mat[,3],PlmPu_pop_mat[,2]-b_len,PlmPu_noise_mat[,3],PlmPu_pop_mat[,2]+b_len,col='dimgrey')
  
  dat<-data.frame(x1=PlmPu_noise_mat[,2],y1=PlmPu_pop_mat[,2])
  mylm<-lm(y1~x1,data=dat)
  abline(mylm,col="black")
  lines(x=xlim1,y=ylim1,col="black",lty="dashed")
  c<-cor.test(dat$x1,dat$y1,method = "pearson",alternative = "t")
  mtext(paste0("Pearson = ",round(unname(c$estimate),3),", p = ",round(c$p.value,3),sep=""),cex=1.5,line=0.1)
  par(op)
  dev.off()
  
  # Plotting D2u of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_D2u_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,D2u_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("D"[u]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.0005+max(D2u_noise_mat[,2],D2u_pop_mat[,2])))
  segments(corcoef_list,D2u_noise_mat[,1],corcoef_list,D2u_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,D2u_noise_mat[,1],corcoef_list+bar_len,D2u_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,D2u_noise_mat[,3],corcoef_list+bar_len,D2u_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,D2u_noise_mat[,1],corcoef_list,D2u_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,D2u_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,D2u_pop_mat[,1],corcoef_list,D2u_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,D2u_pop_mat[,1],corcoef_list+bar_len,D2u_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,D2u_pop_mat[,3],corcoef_list+bar_len,D2u_pop_mat[,3],col='black')
  #arrows(corcoef_list,D2u_pop_mat[,1],corcoef_list,D2u_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_D2u,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting D2l of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_D2l_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,D2l_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("D"[l]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(0,0.0005+max(D2l_noise_mat[,2],D2l_pop_mat[,2])))
  segments(corcoef_list,D2l_noise_mat[,1],corcoef_list,D2l_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,D2l_noise_mat[,1],corcoef_list+bar_len,D2l_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,D2l_noise_mat[,3],corcoef_list+bar_len,D2l_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,D2l_noise_mat[,1],corcoef_list,D2l_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,D2l_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,D2l_pop_mat[,1],corcoef_list,D2l_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,D2l_pop_mat[,1],corcoef_list+bar_len,D2l_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,D2l_pop_mat[,3],corcoef_list+bar_len,D2l_pop_mat[,3],col='black')
  #arrows(corcoef_list,D2l_pop_mat[,1],corcoef_list,D2l_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_D2l,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting D2u - D2l of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_D2u-D2l_vs_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(3.5,4.5,2,3.5), mgp=c(1.9,0.5,0))
  plot(corcoef_list,D2umD2l_noise_mat[,2],cex=0.5,col="dimgrey",xlab=xlabel,ylab=expression("D"["u"]^2-"D"["l"]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=c(0,1),ylim=c(-0.004,0.004)) 
  #ylim=c(0,0.002+max(D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,2])))
  segments(corcoef_list,D2umD2l_noise_mat[,1],corcoef_list,D2umD2l_noise_mat[,3],col='dimgrey')
  segments(corcoef_list-bar_len,D2umD2l_noise_mat[,1],corcoef_list+bar_len,D2umD2l_noise_mat[,1],col='dimgrey')
  segments(corcoef_list-bar_len,D2umD2l_noise_mat[,3],corcoef_list+bar_len,D2umD2l_noise_mat[,3],col='dimgrey')
  #arrows(corcoef_list,D2umD2l_noise_mat[,1],corcoef_list,D2umD2l_noise_mat[,3],length=0.03, angle=90, code=3, col='dimgrey')
  points(corcoef_list,D2umD2l_pop_mat[,2],cex=0.5,col="black")
  segments(corcoef_list,D2umD2l_pop_mat[,1],corcoef_list,D2umD2l_pop_mat[,3],col='black')
  segments(corcoef_list-bar_len,D2umD2l_pop_mat[,1],corcoef_list+bar_len,D2umD2l_pop_mat[,1],col='black')
  segments(corcoef_list-bar_len,D2umD2l_pop_mat[,3],corcoef_list+bar_len,D2umD2l_pop_mat[,3],col='black')
  #arrows(corcoef_list,D2umD2l_pop_mat[,1],corcoef_list,D2umD2l_pop_mat[,3],length=0.03, angle=90, code=3, col='black')
  lines(range(0,1),c(0,0),type='l',lty='dotted',col='black')
  par(new = TRUE)
  plot(corcoef_list,pval_D2umD2l,col="black",type="p",axes = FALSE, bty = "n", xlab = "", ylab = "",
       xlim=c(0,1),ylim=c(0,1),pch=2)
  lines(range(0,1),c(0.05,0.05),type='l',lty='dashed',col='black')
  axis(side=4,col='black',col.axis="black",cex.axis=1.5)
  mtext(side = 4, line = 2, 'p values', col='black',cex=2)
  par(op)
  dev.off()
  
  # Plotting correlation between D2u - D2l of noise and pop copula
  pdf(paste0(resloc,BiCopName(fcode,short=F),"_scatter_D2umD2l_",xlabel,".pdf",sep=""),height=4,width=5)
  op<-par(mar=c(4.5,6,2,2), mgp=c(3,0.5,0))
  ylim1<-range(D2umD2l_pop_mat[,1],D2umD2l_pop_mat[,3])
  xlim1<-range(D2umD2l_noise_mat[,1],D2umD2l_noise_mat[,3])
  
  plot(D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,2],cex=2,col="black",
       xlab=expression("D"["u,input"]^2-"D"["l,input"]^2),
       ylab=expression("D"["u,output"]^2-"D"["l,output"]^2),
       cex.lab=2,cex.axis=1.5,
       xlim=xlim1,ylim=ylim1)
  
  b_len<-diff(xlim1)/50
  segments(D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,1],D2umD2l_noise_mat[,2],D2umD2l_pop_mat[,3],col='black')
  segments(D2umD2l_noise_mat[,2]-b_len,D2umD2l_pop_mat[,1],D2umD2l_noise_mat[,2]+b_len,D2umD2l_pop_mat[,1],col='black')
  segments(D2umD2l_noise_mat[,2]-b_len,D2umD2l_pop_mat[,3],D2umD2l_noise_mat[,2]+b_len,D2umD2l_pop_mat[,3],col='black')
  
  segments(D2umD2l_noise_mat[,1],D2umD2l_pop_mat[,2],D2umD2l_noise_mat[,3],D2umD2l_pop_mat[,2],col='dimgrey')
  segments(D2umD2l_noise_mat[,1],D2umD2l_pop_mat[,2]-b_len,D2umD2l_noise_mat[,1],D2umD2l_pop_mat[,2]+b_len,col='dimgrey')
  segments(D2umD2l_noise_mat[,3],D2umD2l_pop_mat[,2]-b_len,D2umD2l_noise_mat[,3],D2umD2l_pop_mat[,2]+b_len,col='dimgrey')
  
  dat<-data.frame(x1=D2umD2l_noise_mat[,2],y1=D2umD2l_pop_mat[,2])
  mylm<-lm(y1~x1,data=dat)
  abline(mylm,col="black")
  lines(x=xlim1,y=ylim1,col="black",lty="dashed")
  c<-cor.test(dat$x1,dat$y1,method = "pearson",alternative = "t")
  mtext(paste0("Pearson = ",round(unname(c$estimate),3),", p = ",round(c$p.value,3),sep=""),cex=1.5,line=0.1)
  par(op)
  dev.off()
  
  #op2<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  #plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #legend("top", c("noise","population", "p-value(paired t-test)"), col = c("dimgrey", "black", "black"),
  #       cex = 0.8, pch = c(16, 16, 1), xpd = TRUE, horiz = T, inset = c(0,0), 
  #       bty = "n") 
  #par(op2)
  
  pdf(paste0(resloc,"common_legend_cause4copula_stat.pdf",sep=""),height=1,width=15)
  op<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("center", c("input noise copula","model output copula", "p-value(paired t-test)"), col = c("grey70", "black", "black"),
         pt.lwd=4.5, cex = 2.5, pch = c(1, 1, 2), xpd = TRUE, horiz = T, inset = c(0,0), text.col = c("grey70", "black", "black"),
         bty = "n")
  par(op)
  dev.off()
  
  
}