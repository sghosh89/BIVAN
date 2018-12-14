library(VineCopula)
f<-function(par,destau,desLTmUT,fam){
  cop<-BiCop(family=fam,par=par[1],par2=par[2])
  tau<-cop$tau 
  LT<-cop$taildep$lower
  UT<-cop$taildep$upper
  LTmUT<-LT-UT
  z<-(tau-destau)^2
  z<-z+(LTmUT-desLTmUT)^2
  return(z)
}

pdf("./Results/copula_pedagog_fig_BB1.pdf",height=8,width=12)
op<-par(mfrow=c(2,3),mar=c(8,8,0.5,2), mgp=c(5,2,0))

desLTmUT_vector<-c(-0.4,0,0.4)

for(i in c(1:length(desLTmUT_vector))){
  # print(i)
  destau<-0.5
  call_optim<-optim(par=c(0.2,1.5),fn=f,destau=destau,desLTmUT=desLTmUT_vector[i],fam=7,method = "BFGS")
  obj<-BiCop(family=7,par=call_optim$par[1],par2=call_optim$par[2])
  c<-BiCopSim(N=500,obj)
  plot(c[,1],c[,2],col="grey",xlab="u",ylab="v",cex.lab=3.5,cex.axis=3.5,pch=20,cex=1.4)
  text(0.1,0.9,LETTERS[i],cex=6)
  #legend("topleft",LETTERS[i],bty="n",cex=6,x.intersp=-0.5)
  #mtext(paste0("( LT - UT ) = ",round((obj$taildep$lower - obj$taildep$upper),2),sep=""),side=3,cex=2)
  legend("bottomright",c(paste0("KC= ", round(obj$tau,2)),
                         paste0("LT=",round(obj$taildep$lower,2)),
                         paste0("UT=",round(obj$taildep$upper,2))),bty="n",horiz = F,cex=2.5)
  #legend("bottomright",paste0("LT=",round(obj$taildep$lower,2)),bty="n",cex=2)
}

for(i in c(1:length(desLTmUT_vector))){
  # print(i)
  destau<-0.8
  call_optim<-optim(par=c(0.2,1.5),fn=f,destau=destau,desLTmUT=desLTmUT_vector[i],fam=7,method = "BFGS")
  obj<-BiCop(family=7,par=call_optim$par[1],par2=call_optim$par[2])
  c<-BiCopSim(N=500,obj)
  plot(c[,1],c[,2],col="grey",xlab="u",ylab="v",cex.lab=3.5,cex.axis=3.5,pch=20,cex=1.4)
  text(0.1,0.9,LETTERS[3+i],cex=6)
  #mtext(paste0("( LT - UT ) = ",round((obj$taildep$lower - obj$taildep$upper),2),sep=""),side=3,cex=2)
  legend("bottomright",c(paste0("KC= ", round(obj$tau,2)),
                         paste0("LT=",round(obj$taildep$lower,2)),
                         paste0("UT=",round(obj$taildep$upper,2))),bty="n",horiz = F,cex=2.5)
  #legend("bottomright",paste0("LT=",round(obj$taildep$lower,2)),bty="n",cex=2)
}

par(op)
dev.off()