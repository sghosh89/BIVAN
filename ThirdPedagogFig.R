library(copula)
library(mvtnorm)
#set.seed(101)

#***prep the data and copulas

numpts<-500

#this was for reminding myself how the paramaterization of a normal
#copula works
#sig<-matrix(c(1,.7,.7,1),2,2)
#d<-rmvnorm(25000,mean=c(0,0),sigma=sig)
#cor(d[,1],d[,2],method="kendall")
#cor(d[,1],d[,2],method="spearman")

#normal copula corresponding to a covariance of 0.7
nc<-normalCopula(.8,2,dispstr = "un")
dnc<-rCopula(numpts,nc)

#frank copula
fc<-frankCopula(10)
dfc<-rCopula(numpts,fc)

cc<-claytonCopula(5)
sc<-rotCopula(cc)

#***now do the plot

#plot dimensions, units inches
xmarg_ht<-.25
ymarg_wd<-.25
numsp<-.2
gap<-0.1
labgap<-.225
pan_wd<-1
pan_ht<-pan_wd
tot_wd<-ymarg_wd+2*gap+2*pan_wd+numsp
tot_ht<-xmarg_ht+2*gap+2*pan_ht+numsp+2*labgap
pchval<-20
cexvalpts<-0.25
cexvaltxt<-0.75
cexvalpl<-1.5
pdf(file="./Results/PedagogFig3.pdf",width=tot_wd,height=tot_ht)

#upper left, normal pdf
par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+1*pan_ht+1*gap+1*numsp+labgap)/tot_ht,
          (xmarg_ht+2*pan_ht+1*gap+1*numsp+labgap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
x<-seq(from=0,to=1,length.out=101)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)
z<-dCopula(dcarg,nc,log=FALSE)
z<-matrix(z,101,101)
contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),
        xaxt="n")
axis(side=1,labels=FALSE)
mtext("v",2,1.15)
mtext("A",3,.05,cex=cexvalpl,adj=0)

#upper right, normal data
par(fig=c((ymarg_wd+1*numsp+pan_wd+gap)/tot_wd,
          (ymarg_wd+1*numsp+2*pan_wd+gap)/tot_wd,
          (xmarg_ht+1*pan_ht+1*gap+1*numsp+labgap)/tot_ht,
          (xmarg_ht+2*pan_ht+1*gap+1*numsp+labgap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(dnc[,1],dnc[,2],type='p',pch=pchval,cex=cexvalpts,
     xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
axis(side=1,labels=FALSE)
axis(side=2,labels=FALSE)
mtext("B",3,.05,cex=cexvalpl,adj=0)

#lower left, frank pdf
par(fig=c((ymarg_wd+numsp)/tot_wd,
          (ymarg_wd+numsp+pan_wd)/tot_wd,
          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
z<-dCopula(dcarg,fc,log=FALSE)
z<-matrix(z,101,101)
contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),
        xlab="u",ylab="v")
mtext("u",1,1.15)
mtext("v",2,1.15)
mtext("C",3,.05,cex=cexvalpl,adj=0)

#lower right, frank data
par(fig=c((ymarg_wd+1*numsp+pan_wd+gap)/tot_wd,
          (ymarg_wd+1*numsp+2*pan_wd+gap)/tot_wd,
          (xmarg_ht+0*pan_ht+0*gap+1*numsp)/tot_ht,
          (xmarg_ht+1*pan_ht+0*gap+1*numsp)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(dfc[,1],dfc[,2],type='p',pch=pchval,cex=cexvalpts,
     xlim=c(0,1),ylim=c(0,1),yaxt="n")
axis(side=2,labels=FALSE)
mtext("u",1,1.15)
mtext("D",3,.05,cex=cexvalpl,adj=0)

dev.off()