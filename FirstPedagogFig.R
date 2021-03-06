library(copula)

#***prep the data

#clayton data
cc5<-claytonCopula(5)
dcc5<-rCopula(250,cc5)
#plot(dcc5[,1],dcc5[,2],type='p')

#survival clayton data
sc5<-rotCopula(cc5)
#dsc5<-rCopula(250,sc5)
dsc5<-(-dcc5+1)
#plot(dsc5[,1],dsc5[,2],type='p')

#make marginals normal
dcc5_nm<-qnorm(dcc5)
#plot(dcc5_nm[,1],dcc5_nm[,2],type='p')

dsc5_nm<-qnorm(dsc5)
#plot(dsc5_nm[,1],dsc5_nm[,2],type='p')

#make marginals gamma
dcc5_gm<-qgamma(dcc5,shape=2,scale=2)
#plot(dcc5_gm[,1],dcc5_gm[,2],type='p')

dsc5_gm<-qgamma(dsc5,shape=2,scale=2)
#plot(dsc5_gm[,1],dsc5_gm[,2],type='p')

#normal copula 
nc<-normalCopula(.8,2,dispstr = "un")
ncpar<-iRho(nc,rho(cc5))
nc<-normalCopula(ncpar,2,dispstr = "un")
dnc<-rCopula(250,nc)

#***now do the plot

#plot dimensions, units inches
xmarg_ht<-.25
ymarg_wd<-.25
numsp<-.2
gap<-0.1
pan_wd<-1
pan_ht<-pan_wd
tot_wd<-ymarg_wd+2*gap+2*pan_wd+2*numsp
tot_ht<-xmarg_ht+3*gap+3*pan_ht+3*numsp
pchval<-20
cexvalpts<-0.25
cexvaltxt<-0.75
cexvalpl<-1.5
cexaxs<-0.9
pdf(file="./Results/PedagogFig1.pdf",width=tot_wd,height=tot_ht)

#upper-left panel, normal marginals, left-tail dependence
par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+2*pan_ht+2*gap+3*numsp)/tot_ht,
          (xmarg_ht+3*pan_ht+2*gap+3*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
rg<-max(abs(dcc5_nm))
plot(dcc5_nm[,1],dcc5_nm[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(-rg,rg),ylim=c(-rg,rg),col="grey",cex.axis=cexaxs)
mtext(expression(z[v]),2,line=1.15)
text(-rg,rg,"A",cex=cexvalpl,adj=c(0,1))
text(rg,-rg+.25*rg,paste0("P=",round(cor(dcc5_nm[,1],dcc5_nm[,2]),2)),
     adj=c(1,0),cex=cexvaltxt)
text(rg,-rg,paste0("S=",round(cor(dcc5_nm[,1],dcc5_nm[,2],method="spearman"),2)),
     adj=c(1,0),cex=cexvaltxt)

#upper-right panel, normal marginals, right-tail dependence
par(fig=c((ymarg_wd+2*numsp+pan_wd+gap)/tot_wd,
          (ymarg_wd+2*numsp+2*pan_wd+gap)/tot_wd,
          (xmarg_ht+2*pan_ht+2*gap+3*numsp)/tot_ht,
          (xmarg_ht+3*pan_ht+2*gap+3*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
rg<-max(abs(dsc5_nm))
plot(dsc5_nm[,1],dsc5_nm[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(-rg,rg),ylim=c(-rg,rg),col="grey",cex.axis=cexaxs)
text(-rg,rg,"B",cex=cexvalpl,adj=c(0,1))
text(rg,-rg+.25*rg,paste0("P=",round(cor(dsc5_nm[,1],dsc5_nm[,2]),2)),
     adj=c(1,0),cex=cexvaltxt)
text(rg,-rg,paste0("S=",round(cor(dsc5_nm[,1],dsc5_nm[,2],method="spearman"),2)),
     adj=c(1,0),cex=cexvaltxt)

#2nd row of panels, left, uniform marginals, left-tail dependence
par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+1*pan_ht+1*gap+2*numsp)/tot_ht,
          (xmarg_ht+2*pan_ht+1*gap+2*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
#the data
rg1<-0
rg2<-1
dcc5pobs<-copula::pobs(dcc5)
plot(dcc5pobs[,1],dcc5pobs[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey",cex.axis=cexaxs)
mtext(expression(z[v]),2,line=1.15)
#the contours
#x<-seq(from=0,to=1,length.out=101)
#y<-x
#dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)
#z<-dCopula(dcarg,cc5)
#z<-matrix(z,101,101)
#contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),col="grey",labels ="",add=T)
#contour(x,y,log10(z),nlevels=5,lty=0, add=T)
#contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1))
text(rg1,rg2,"C",cex=cexvalpl,adj=c(0,1))
text(rg2,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dcc5pobs[,1],dcc5pobs[,2]),2)),
     adj=c(1,0),cex=cexvaltxt)
text(rg2,rg1,paste0("S=",round(cor(dcc5pobs[,1],dcc5pobs[,2],method="spearman"),2)),
     adj=c(1,0),cex=cexvaltxt)

#2nd row of panels, right, uniform marginals, right-tail dependence
par(fig=c((ymarg_wd+2*numsp+pan_wd+gap)/tot_wd,
          (ymarg_wd+2*numsp+2*pan_wd+gap)/tot_wd,
          (xmarg_ht+1*pan_ht+1*gap+2*numsp)/tot_ht,
          (xmarg_ht+2*pan_ht+1*gap+2*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
#the data
dsc5pobs<-copula::pobs(dsc5)
plot(dsc5pobs[,1],dsc5pobs[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey",cex.axis=cexaxs)
#the contours
#x<-seq(from=0,to=1,length.out=101)
#y<-x
#dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)
#z<-dCopula(dcarg,sc5)
#z<-matrix(z,101,101)
#contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),col="grey",labels ="",add=T)
#contour(x,y,log10(z),nlevels=5,lty=0, add=T)
text(rg1,rg2,"D",cex=cexvalpl,adj=c(0,1))
text(rg2,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dsc5pobs[,1],dsc5pobs[,2]),2)),
     adj=c(1,0),cex=cexvaltxt)
text(rg2,rg1,paste0("S=",round(cor(dsc5pobs[,1],dsc5pobs[,2],method="spearman"),2)),
     adj=c(1,0),cex=cexvaltxt)

#3th row of panels, left
par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
dcc5fv<-dcc5pobs
dcc5fv[,1]<-(-dcc5fv[,1]+1)
rg1<-0
rg2<-1
plot(dcc5fv[,1],dcc5fv[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey",cex.axis=cexaxs)
text(rg2,rg2,"E",cex=cexvalpl,adj=c(1,1))
text(rg1,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dcc5fv[,1],dcc5fv[,2]),2)),
     adj=c(0,0),cex=cexvaltxt)
text(rg1,rg1,paste0("S=",round(cor(dcc5fv[,1],dcc5fv[,2],method="spearman"),2)),
     adj=c(0,0),cex=cexvaltxt)
mtext(expression(z[v]),2,line=1.15)
mtext(expression(z[h]),1,line=1.15)

#3th row of panels, right
par(fig=c((ymarg_wd+2*numsp+pan_wd+gap)/tot_wd,
          (ymarg_wd+2*numsp+2*pan_wd+gap)/tot_wd,
          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
dsc5fv<-dsc5pobs
dsc5fv[,1]<-(-dsc5fv[,1]+1)
plot(dsc5fv[,1],dsc5fv[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey",cex.axis=cexaxs)
text(rg2,rg2,"F",cex=cexvalpl,adj=c(1,1))
text(rg1,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dsc5fv[,1],dsc5fv[,2]),2)),
     adj=c(0,0),cex=cexvaltxt)
text(rg1,rg1,paste0("S=",round(cor(dsc5fv[,1],dsc5fv[,2],method="spearman"),2)),
     adj=c(0,0),cex=cexvaltxt)
mtext(expression(z[h]),1,line=1.15)

#4th row of panels, left - normal copula data and contours
#par(fig=c((ymarg_wd+numsp)/tot_wd,
#          (ymarg_wd+numsp+pan_wd)/tot_wd,
#          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
#          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
#    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
#data
#plot(dnc[,1],dnc[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey")
##text(rg2,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dnc[,1],dnc[,2]),2)),
##     adj=c(1,0),cex=cexvaltxt)
##text(rg2,rg1,paste0("S=",round(cor(dnc[,1],dnc[,2],method="spearman"),2)),
##     adj=c(1,0),cex=cexvaltxt)
##contours
#x<-seq(from=0,to=1,length.out=101)
#y<-x
#dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)
#z<-dCopula(dcarg,nc,log=FALSE)
#z<-matrix(z,101,101)
#contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),xaxt="n",col="grey",labels ="",add=T)
##contour(x,y,log10(z),nlevels=5,lty=0, add=T)
#axis(side=1,labels=FALSE)
#text(rg1,rg2,"G",cex=cexvalpl,adj=c(0,1))
#mtext(expression(z[h]),1,line=1.15)
#mtext(expression(z[v]),2,line=1.15)

##4th row of panels, right - figure to get cor_lb,ub across
#par(fig=c((ymarg_wd+2*numsp+pan_wd+gap)/tot_wd,
#          (ymarg_wd+2*numsp+2*pan_wd+gap)/tot_wd,
#          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
#          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
#    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
#plot(dcc5[,1],dcc5[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey")
#lines(c(0,.4),c(.4,0),type='l')
#lines(c(0,.65),c(.65,0),type='l')
#lines(c(0,1.6),c(1.6,0),type='l')
#lines(c(0,1.35),c(1.35,0),type='l')
#mtext(expression(z[h]),1,line=1.15)
#text(rg1,rg2,"H",cex=cexvalpl,adj=c(0,1))

dev.off()

#***A pedagog fig for intro-ing the nonparam stats

tot_wd<-ymarg_wd+numsp+pan_wd+gap
tot_ht<-xmarg_ht+numsp+pan_ht+gap
pdf(file="./Results/NonparamStatsFig.pdf",width=tot_wd,height=tot_ht)

par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+numsp)/tot_ht,
          (xmarg_ht+numsp+pan_ht)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
plot(dcc5[,1],dcc5[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey",cex.axis=cexaxs,xaxs="i",yaxs="i")
lines(c(0,.4),c(.4,0),type='l')
lines(c(0,.65),c(.65,0),type='l')
lines(c(0,1.6),c(1.6,0),type='l')
lines(c(0,1.35),c(1.35,0),type='l')
mtext("u",1,line=1.15)
mtext("v",2,line=1.15)
#text(rg1,rg2,"H",cex=cexvalpl,adj=c(0,1))

dev.off()

#***a pedagog fig for the sup mat about monotonic transformations

#plot dimensions, units inches
xmarg_ht<-.25
ymarg_wd<-.25
numsp<-.2
gap<-0.1
pan_wd<-2
pan_ht<-pan_wd
tot_wd<-ymarg_wd+2*gap+2*pan_wd+2*numsp
tot_ht<-xmarg_ht+1*gap+1*pan_ht+1*numsp
pchval<-20
cexvalpts<-0.25
cexvaltxt<-0.75
cexvalpl<-1.5
pdf(file="./Results/PedagogFig_SMExport.pdf",width=tot_wd,height=tot_ht)

#panel left, gamma marginals, left-tail dependence
par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
rg1<-min(dcc5_gm,dsc5_gm)
rg2<-max(dcc5_gm,dsc5_gm)
plot(dcc5_gm[,1],dcc5_gm[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey")
mtext(expression(z[v]),2,line=1.15)
text(rg1,rg2,"A",cex=cexvalpl,adj=c(0,1))
text(rg2,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dcc5_gm[,1],dcc5_gm[,2]),2)),
     adj=c(1,0),cex=cexvaltxt)
text(rg2,rg1,paste0("S=",round(cor(dcc5_gm[,1],dcc5_gm[,2],method="spearman"),2)),
     adj=c(1,0),cex=cexvaltxt)
mtext(expression(z[h]),1,line=1.15)

#panel right, gamma marginals, right-tail dependence
par(fig=c((ymarg_wd+2*numsp+pan_wd+gap)/tot_wd,
          (ymarg_wd+2*numsp+2*pan_wd+gap)/tot_wd,
          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=T)
plot(dsc5_gm[,1],dsc5_gm[,2],type='p',pch=pchval,cex=cexvalpts,xlim=c(rg1,rg2),ylim=c(rg1,rg2),col="grey")
text(rg1,rg2,"B",cex=cexvalpl,adj=c(0,1))
text(rg2,rg1+.125*(rg2-rg1),paste0("P=",round(cor(dsc5_gm[,1],dsc5_gm[,2]),2)),
     adj=c(1,0),cex=cexvaltxt)
text(rg2,rg1,paste0("S=",round(cor(dsc5_gm[,1],dsc5_gm[,2],method="spearman"),2)),
     adj=c(1,0),cex=cexvaltxt)
mtext(expression(z[h]),1,line=1.15)

dev.off()


