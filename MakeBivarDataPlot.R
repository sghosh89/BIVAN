#***plotting dimensions, units inches
xmarg_ht<-.45
ymarg_wd<-.45
gap<-0.1
tot_wd<-6.5
pan_wd<-(tot_wd-4*ymarg_wd-4*gap)/4
pan_ht<-pan_wd
tot_ht<-2*xmarg_ht+2*pan_ht+2*gap
pchval<-20
cexvalpts<-0.25
cexvaltxt<-0.75
cexvalpl<-1.5
pdf(file="./Results/BivarDataPlot.pdf",width=tot_wd,height=tot_ht)

#***top row of panels

#panel 1
px<-0
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
#Shya, please load the appropriate data here as d1, replacing the below 
#nonsense filler row
d1<-data.frame(x=1:3,y=1:3)
d1c<-copula::pobs(d1)
plot(d1$x,d1$y,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("Soil C",1,1.15)
mtext("Soil N",2,1.15)
mtext("A",3,-1.4,cex=cexvalpl,adj=.05)

#panel 2
px<-1
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
#Shya, please load the appropriate data here as d2, replacing the below 
#nonsense filler row. Please make sure you use log10 here.
d2<-data.frame(x=1:3,y=1:3)
d2c<-copula::pobs(d2)
plot(d2$x,d2$y,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("log(body mass)",1,1.15)
mtext("log(BMR)",2,1.15)
mtext("B",3,-1.4,cex=cexvalpl,adj=.05)

#panel 3
px<-2
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
#Shya, please load the appropriate data here as d3, replacing the below 
#nonsense filler row. Please make sure you use log10 here.
d3<-data.frame(x=1:3,y=1:3)
d3c<-copula::pobs(d3)
plot(d3$x,d3$y,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("log(body mass)",1,1.15)
mtext("log(BMR)",2,1.15)
mtext("C",3,-1.4,cex=cexvalpl,adj=.05)

#panel 4
px<-3
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
#Shya, please load the appropriate data here as d4, replacing the below 
#nonsense filler row. 
d4<-data.frame(x=1:3,y=1:3)
d4c<-copula::pobs(d4)
plot(d4$x,d4$y,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("Avg. H",1,1.15)
mtext("Avg. biomass",2,1.15)
mtext("D",3,-1.4,cex=cexvalpl,adj=.05)

#***bottom row of panels

#panel 1
px<-0
py<-0
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d1c[,1],d1c[,2],type='p',pch=pchval,cex=cexvalpts,
     col="grey",xlab="u",ylab="v",xlim=c(0,1),ylim=c(0,1))
mtext("u",1,1.15)
mtext("v",2,1.15)
mtext("E",3,-1.4,cex=cexvalpl,adj=.05)

#panel 2
px<-1
py<-0
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d2c[,1],d2c[,2],type='p',pch=pchval,cex=cexvalpts,
     col="grey",xlab="u",ylab="v",xlim=c(0,1),ylim=c(0,1))
mtext("u",1,1.15)
mtext("v",2,1.15)
mtext("F",3,-1.4,cex=cexvalpl,adj=.05)

#panel 3
px<-2
py<-0
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d3c[,1],d3c[,2],type='p',pch=pchval,cex=cexvalpts,
     col="grey",xlab="u",ylab="v",xlim=c(0,1),ylim=c(0,1))
mtext("u",1,1.15)
mtext("v",2,1.15)
mtext("G",3,-1.4,cex=cexvalpl,adj=.05)

#panel 4
px<-3
py<-0
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(d4c[,1],d4c[,2],type='p',pch=pchval,cex=cexvalpts,
     col="grey",xlab="u",ylab="v",xlim=c(0,1),ylim=c(0,1))
mtext("u",1,1.15)
mtext("v",2,1.15)
mtext("H",3,-1.4,cex=cexvalpl,adj=.05)

dev.off()
