#***
#This fle generates a pedagogical figure for showing people how different copulas
#can be combined with different marginals to give different distributions
#***

library(copula)

#***make the two copulas

#**Gaussian copula
ncop<-normalCopula(.7)

#**Clayton copula with same Spearman
ccop<-claytonCopula(2)
ccop<-claytonCopula(iRho(ccop,rho(ncop)))

#***make the marginals
dnmarg<-function(x){return(dnorm(x,mean=0,sd=1))}
pnmarg<-function(x){return(pnorm(x,mean=0,sd=1))}
qnmarg<-function(x){return(qnorm(x,mean=0,sd=1))}

shape<-5
scale<-1
dgmarg<-function(x){return(dgamma(x,shape=shape,scale=scale))}
pgmarg<-function(x){return(pgamma(x,shape=shape,scale=scale))}
qgmarg<-function(x){return(qgamma(x,shape=shape,scale=scale))}

#***other stuff
numsamps<-50

#***make the plot

#**set the plot up, units inches
totwd<-6.5
biggap<-.2
smallgap<-.1
xmarht<-.5
ymarwd<-xmarht
textspace<-0.2
panwd<-(totwd-3*ymarwd-2*biggap-4*smallgap-textspace)/(4)
margpansm<-panwd/3
panht<-panwd
margpanbg<-panht
totht<-2*xmarht+2*panht+2*margpansm+biggap+3*smallgap+textspace
if (!isTRUE(all.equal(totwd,3*ymarwd+3*panwd+3*margpansm+2*biggap+4*smallgap+textspace)))
{
  stop("Error in figure dimensions")
}
pdf(file="./Results/PedagogFigCopThry.pdf",width=totwd,height=totht)

#**panel for contours for the normal cop
par(fig=c((textspace+ymarwd)/totwd,
          (textspace+ymarwd+panwd)/totwd,
          (2*xmarht+panht+margpansm+smallgap+biggap)/totht,
          (2*xmarht+2*panht+margpansm+smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)

x<-seq(from=0,to=1,length.out=101)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)

z<-dCopula(dcarg,ncop,log=FALSE)
z<-matrix(z,101,101)
levs<-c(-3,-2,-1,0)
contour(x,y,log10(z),levels=levs,xlim=c(0,1),ylim=c(0,1),col="gray",labels="")

mtext("u",1,1.2)
mtext("v",2,1.2)
mtext("A",3,.8,at=c(-.35))
mtext("NORMAL COPULA",2,2.4)

#***add samples
set.seed(101)
samps_ncop<-rCopula(numsamps,ncop)
points(samps_ncop[,1],samps_ncop[,2],type="p",pch=20,col="grey",cex=.5)

contour(x,y,log10(z),levels=levs,lty=0,add=TRUE)

#**normal copula plot, top marginal
par(fig=c((textspace+ymarwd)/totwd,
          (textspace+ymarwd+margpanbg)/totwd,
          (2*xmarht+2*panht+margpansm+2*smallgap+biggap)/totht,
          (2*xmarht+2*panht+2*margpansm+2*smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=1,by=0.01)
y<-rep(1,length(x))

plot(x,y,type='l',yaxt="n",xaxt="n")
rug(samps_ncop[,1],ticksize=0.06)
axis(side=1,labels=FALSE)
mtext("UNIFORM MARGINALS",3,0.2)

#**normal copula plot, right marginal
par(fig=c((textspace+ymarwd+panwd+smallgap)/totwd,
          (textspace+ymarwd+panwd+smallgap+margpansm)/totwd,
          (2*xmarht+panht+margpansm+smallgap+biggap)/totht,
          (2*xmarht+2*panht+margpansm+smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

y<-seq(from=0,to=1,by=0.01)
x<-rep(1,length(y))

plot(x,y,type='l',xaxt="n",yaxt="n")
rug(samps_ncop[,2],ticksize=0.06,side=2)
axis(side=2,labels=FALSE)

#**first normal copula example, main panel
par(fig=c((textspace+2*ymarwd+panwd+biggap+margpansm)/totwd,
          (textspace+2*ymarwd+2*panwd+biggap+margpansm)/totwd,
          (2*xmarht+panht+margpansm+smallgap+biggap)/totht,
          (2*xmarht+2*panht+margpansm+smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=-3,to=3,by=0.05)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)

z<-dCopula(cbind(pnmarg(dcarg[,1]),pnmarg(dcarg[,2])),ncop,log=FALSE)*dnmarg(dcarg[,1])*dnmarg(dcarg[,2])
z<-matrix(z,length(x),length(y))
levs<-c(-8,-6,-4,-2,-1)
contour(x,y,log10(z),xlim=c(-3,3),ylim=c(-3,3),levels=levs,col="gray",labels="")

mtext("x",1,1.2)
mtext("y",2,1.2)
mtext("C",3,.8,at=c(-3-6*.35))

#**add samples
samps<-qnmarg(samps_ncop)
points(samps[,1],samps[,2],type='p',pch=20,cex=0.5,col="grey")

contour(x,y,log10(z),levels=levs,lty=0,add=TRUE)

#**first normal copula example, top marginal
par(fig=c((textspace+2*ymarwd+panwd+biggap+margpansm)/totwd,
          (textspace+2*ymarwd+2*panwd+biggap+margpansm)/totwd,
          (2*xmarht+2*panht+margpansm+2*smallgap+biggap)/totht,
          (2*xmarht+2*panht+2*margpansm+2*smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=-3,to=3,by=0.01)
y<-dnmarg(x)

plot(x,y,type='l',yaxt="n",xaxt="n")
rug(samps[,1],ticksize=0.06)
axis(side=1,labels=FALSE)
mtext("NORMAL MARGINALS",3,0.2)

#**first normal copula example, right marginal
par(fig=c((textspace+2*ymarwd+2*panwd+biggap+smallgap+margpansm)/totwd,
          (textspace+2*ymarwd+2*panwd+biggap+smallgap+margpansm+margpansm)/totwd,
          (2*xmarht+panht+margpansm+smallgap+biggap)/totht,
          (2*xmarht+2*panht+margpansm+smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

y<-seq(from=-3,to=3,by=0.01)
x<-dnmarg(y)

plot(x,y,type='l',xaxt="n",yaxt="n")
rug(samps[,2],ticksize=0.06,side=2)
axis(side=2,labels=FALSE)

#**second normal copula example, main panel
par(fig=c((textspace+3*ymarwd+2*panwd+2*biggap+smallgap+margpansm+margpansm)/totwd,
          (textspace+3*ymarwd+3*panwd+2*biggap+smallgap+margpansm+margpansm)/totwd,
          (2*xmarht+panht+margpansm+smallgap+biggap)/totht,
          (2*xmarht+2*panht+margpansm+smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=16,by=0.05)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)

z<-dCopula(cbind(pgmarg(dcarg[,1]),pgmarg(dcarg[,2])),ncop,log=FALSE)*dgmarg(dcarg[,1])*dgmarg(dcarg[,2])
z<-matrix(z,length(x),length(y))
levs<-c(-10,-5,-4,-3,-2)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,col="gray",labels="")

mtext("x",1,1.2)
mtext("y",2,1.2)
mtext("E",3,.8,at=c(-3-6*.35))

#**add samples
samps<-qgmarg(samps_ncop)
points(samps[,1],samps[,2],type='p',pch=20,cex=0.5,col="grey")

contour(x,y,log10(z),lty=0,add=TRUE,levels=levs)

#**second normal copula example, top marginal
par(fig=c((textspace+3*ymarwd+2*panwd+2*biggap+smallgap+margpansm+margpansm)/totwd,
          (textspace+3*ymarwd+2*panwd+2*biggap+smallgap+margpansm+margpanbg+margpansm)/totwd,
          (2*xmarht+2*panht+margpansm+2*smallgap+biggap)/totht,
          (2*xmarht+2*panht+2*margpansm+2*smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=16,by=0.05)
y<-dgmarg(x)

plot(x,y,type='l',yaxt="n",xaxt="n")
rug(samps[,1],ticksize=0.06)
axis(side=1,labels=FALSE)
mtext("GAMMA MARGINALS",3,0.2)

#**second normal copula example, right marginal
par(fig=c((textspace+3*ymarwd+3*panwd+2*biggap+2*smallgap+margpansm+margpansm)/totwd,
          (textspace+3*ymarwd+3*panwd+2*biggap+2*smallgap+2*margpansm+margpansm)/totwd,
          (2*xmarht+panht+margpansm+smallgap+biggap)/totht,
          (2*xmarht+2*panht+margpansm+smallgap+biggap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

y<-seq(from=0,to=16,by=0.05)
x<-dgmarg(y)

plot(x,y,type='l',xaxt="n",yaxt="n")
rug(samps[,2],ticksize=0.06,side=2)
axis(side=2,labels=FALSE)

#**panel for contours for the clayton cop
par(fig=c((textspace+ymarwd)/totwd,
          (textspace+ymarwd+panwd)/totwd,
          (xmarht)/totht,
          (xmarht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=1,length.out=101)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)

z<-dCopula(dcarg,ccop,log=FALSE)
z<-matrix(z,101,101)
levs<-c(-1.5,-1,-.5,0,.5)
contour(x,y,log10(z),levels=levs,xlim=c(0,1),ylim=c(0,1),col="gray",labels="")

mtext("u",1,1.2)
mtext("v",2,1.2)
mtext("B",3,.8,at=c(-.35))
mtext("CLAYTON COPULA",2,2.4)

#***add samples
set.seed(107)
samps_ccop<-rCopula(numsamps,ccop)
points(samps_ccop[,1],samps_ccop[,2],type="p",pch=20,col="grey",cex=.5)

contour(x,y,log10(z),levels=levs,lty=0,add=TRUE)

#**clayton copula plot, top marginal
par(fig=c((textspace+ymarwd)/totwd,
          (textspace+ymarwd+margpanbg)/totwd,
          (xmarht+panht+smallgap)/totht,
          (xmarht+panht+smallgap+margpansm)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=1,by=0.01)
y<-rep(1,length(x))

plot(x,y,type='l',yaxt="n",xaxt="n")
rug(samps_ccop[,1],ticksize=0.06)
axis(side=1,labels=FALSE)

#**clayton copula plot, right marginal
par(fig=c((textspace+ymarwd+panwd+smallgap)/totwd,
          (textspace+ymarwd+panwd+smallgap+margpansm)/totwd,
          (xmarht)/totht,
          (xmarht+margpanbg)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

y<-seq(from=0,to=1,by=0.01)
x<-rep(1,length(y))

plot(x,y,type='l',xaxt="n",yaxt="n")
rug(samps_ccop[,2],ticksize=0.06,side=2)
axis(side=2,labels=FALSE)

#**first clayton copula example, main panel
par(fig=c((textspace+2*ymarwd+panwd+biggap+margpansm)/totwd,
          (textspace+2*ymarwd+2*panwd+biggap+margpansm)/totwd,
          (1*xmarht)/totht,
          (1*xmarht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=-3,to=3,by=0.05)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)

z<-dCopula(cbind(pnmarg(dcarg[,1]),pnmarg(dcarg[,2])),ccop,log=FALSE)*dnmarg(dcarg[,1])*dnmarg(dcarg[,2])
z<-matrix(z,length(x),length(y))
contour(x,y,log10(z),xlim=c(-3,3),ylim=c(-3,3),nlevels=8,col="gray",labels="")

mtext("x",1,1.2)
mtext("y",2,1.2)
mtext("D",3,.8,at=c(-3-6*.35))

#**add samples
samps<-qnmarg(samps_ccop)
points(samps[,1],samps[,2],type='p',pch=20,cex=0.5,col="grey")

contour(x,y,log10(z),nlevels=8,lty=0,add=TRUE)

#**first clayton copula example, top marginal
par(fig=c((textspace+2*ymarwd+panwd+biggap+margpansm)/totwd,
          (textspace+2*ymarwd+2*panwd+biggap+margpansm)/totwd,
          (xmarht+panht+smallgap)/totht,
          (xmarht+panht+smallgap+margpansm)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=-3,to=3,by=0.01)
y<-dnmarg(x)

plot(x,y,type='l',yaxt="n",xaxt="n")
rug(samps[,1],ticksize=0.06)
axis(side=1,labels=FALSE)

#**first clayton copula example, right marginal
par(fig=c((textspace+2*ymarwd+2*panwd+biggap+smallgap+margpansm)/totwd,
          (textspace+2*ymarwd+2*panwd+biggap+smallgap+margpansm+margpansm)/totwd,
          (xmarht)/totht,
          (xmarht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

y<-seq(from=-3,to=3,by=0.01)
x<-dnmarg(y)

plot(x,y,type='l',xaxt="n",yaxt="n")
rug(samps[,2],ticksize=0.06,side=2)
axis(side=2,labels=FALSE)

#**second clayton copula example, main panel
par(fig=c((textspace+3*ymarwd+2*panwd+2*biggap+smallgap+margpansm+margpansm)/totwd,
          (textspace+3*ymarwd+3*panwd+2*biggap+smallgap+margpansm+margpansm)/totwd,
          (xmarht)/totht,
          (xmarht+panht)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=16,by=0.05)
y<-x
dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)

z<-dCopula(cbind(pgmarg(dcarg[,1]),pgmarg(dcarg[,2])),ccop,log=FALSE)*dgmarg(dcarg[,1])*dgmarg(dcarg[,2])
z<-matrix(z,length(x),length(y))
levs<-c(-10,-5,-4,-3,-2,-1)
contour(x,y,log10(z),xlim=range(x),ylim=range(y),levels=levs,col="gray",labels="")

mtext("x",1,1.2)
mtext("y",2,1.2)
mtext("F",3,.8,at=c(-3-6*.35))

#**add samples
samps<-qgmarg(samps_ccop)
points(samps[,1],samps[,2],type='p',pch=20,cex=0.5,col="grey")

contour(x,y,log10(z),lty=0,add=TRUE,levels=levs)

#**second clayton copula example, top marginal
par(fig=c((textspace+3*ymarwd+2*panwd+2*biggap+smallgap+margpansm+margpansm)/totwd,
          (textspace+3*ymarwd+2*panwd+2*biggap+smallgap+margpansm+margpanbg+margpansm)/totwd,
          (xmarht+panht+smallgap)/totht,
          (xmarht+panht+margpansm+smallgap)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

x<-seq(from=0,to=16,by=0.05)
y<-dgmarg(x)

plot(x,y,type='l',yaxt="n",xaxt="n")
rug(samps[,1],ticksize=0.06)
axis(side=1,labels=FALSE)

#**second clayton copula example, right marginal
par(fig=c((textspace+3*ymarwd+3*panwd+2*biggap+2*smallgap+margpansm+margpansm)/totwd,
          (textspace+3*ymarwd+3*panwd+2*biggap+2*smallgap+2*margpansm+margpansm)/totwd,
          (xmarht)/totht,
          (xmarht+margpanbg)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)

y<-seq(from=0,to=16,by=0.05)
x<-dgmarg(y)

plot(x,y,type='l',xaxt="n",yaxt="n")
rug(samps[,2],ticksize=0.06,side=2)
axis(side=2,labels=FALSE)

dev.off()
