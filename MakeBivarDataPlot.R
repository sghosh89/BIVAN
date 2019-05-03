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

#panel 1 : raw data plot for SoilCN
px<-0
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
d1<-readRDS("./Data/dataForShymolinaCoord.rds")
d1<-d1[,c("SOCstock100","TSNstock100")]
d1<-na.omit(d1)
d1c<-copula::pobs(d1)
plot(log10(d1$SOCstock100),log10(d1$TSNstock100),type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("log(Soil C)",1,1.15)
mtext("log(Soil N)",2,1.15)
mtext("A",3,-1.4,cex=cexvalpl,adj=.05)

#panel 2 : raw data plot for Birds' BMR
px<-1
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
raw_b<-readRDS("./Data/BMR/BirdBodyMassesMetabolicRates/my_birdsBMR_data.RDS")
d2<-log10(raw_b)
colnames(d2)<-c("log10M(g)","log10BMR(KJ/h)")
d2c<-copula::pobs(d2)
plot(d2$`log10M(g)`,d2$`log10BMR(KJ/h)`,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("log(body mass)",1,1.15)
mtext("log(BMR)",2,1.15)
mtext("B",3,-1.4,cex=cexvalpl,adj=.05)

#panel 3 : raw data plot for mammals' BMR
px<-2
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
raw_m<-readRDS("./Data/BMR/MammalianBodyMassesMetabolicRates/my_mammalsBMR_data.RDS")
d3<-log10(raw_m)
colnames(d3)<-c("log10M(g)","log10BMR(KJ/h)")
d3<-as.data.frame(d3)
d3c<-copula::pobs(d3)
plot(d3$`log10M(g)`,d3$`log10BMR(KJ/h)`,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("log(body mass)",1,1.15)
mtext("log(BMR)",2,1.15)
mtext("C",3,-1.4,cex=cexvalpl,adj=.05)

#panel 4 : raw data plot for cedercreek_2000
px<-3
py<-1
par(fig=c(((px+1)*ymarg_wd+px*pan_wd+px*gap)/tot_wd,
          ((px+1)*ymarg_wd+(px+1)*pan_wd+px*gap)/tot_wd,
          ((py+1)*xmarg_ht+py*pan_ht+py*gap)/tot_ht,
          ((py+1)*xmarg_ht+(py+1)*pan_ht+py*gap)/tot_ht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
sr<-readRDS("./Results/Ceder_creek_results/raw_cop_cc_2000.RDS")
d4<-sr$tab_dbp[,c(2:3)]
d4c<-sr$cop_tab_dbp[,c(2:3)]
plot(d4$avg.H,d4$avg.Biomass,type='p',pch=pchval,cex=cexvalpts,
     col="grey")
mtext("H",1,1.15)
mtext("Biomass",2,1.15)
mtext("D",3,-1.4,cex=cexvalpl,adj=.05)

#***bottom row of panels : copula plots of corresponding figs of upper panel

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
