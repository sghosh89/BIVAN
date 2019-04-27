library(VineCopula)
source("CopulaFunctions_flexible.R")
#====================================================
# function to generate copula plot and save results
# Args :
#     nametag : a character variable 
#     dataloc : the location folder to read the data
#     resloc : the location folder to save the results
#
# Output :
#        A plot with raw and copula
#        A list with 100 copula

get_raw_cop<-function(nametag,dataloc,resloc){
  
  pdf(paste(resloc,"raw_and_copula_plot_",nametag,".pdf",sep=""),height=45,width=45)
  op<-par(mfrow=c(15,15))
  
  
  cop_list<-vector("list",100)
  for(i in 1:100){
    d<-read.csv(paste(dataloc,"rep",i-1,".csv",sep=""))
    d<-d[,-1]
    
    plot(d[,1],d[,2],col=rgb(0,0,0,0.2),pch=19,xlab=paste("raw_rep_",i-1,"_T1",sep=""),ylab=paste("raw_rep_",i-1,"_T2",sep=""))
    
    v<-VineCopula::pobs(d)
    cop_list[[i]]<-v
    plot(v[,1],v[,2],col=rgb(1,0,0,0.2),pch=19,xlab=paste("cop_rep_",i-1,"_T1",sep=""),ylab=paste("cop_rep_",i-1,"_T2",sep=""))
  }
  
  par(op)
  dev.off()
  
  return(cop_list)
}


#----------call the function to save copula results---------

if(!dir.exists("./Results/BMR_results/evo_cop")){
  dir.create("./Results/BMR_results/evo_cop")
}
resloc<-"./Results/BMR_results/evo_cop/"

dataloc<-"./Data/BMR/simulations/sims/BivariateNormal/"
nametag<-"BivNormal"
ans<-get_raw_cop(nametag=nametag,dataloc=dataloc,resloc=resloc)
saveRDS(ans,paste(resloc,nametag,"_coplist.RDS",sep=""))

dataloc<-"./Data/BMR/simulations/sims/ExtremeLeftTailDep/"
nametag<-"ExtremeL"
ans<-get_raw_cop(nametag=nametag,dataloc=dataloc,resloc=resloc)
saveRDS(ans,paste(resloc,nametag,"_coplist.RDS",sep=""))

dataloc<-"./Data/BMR/simulations/sims/ExtremeRightTailDep/"
nametag<-"ExtremeR"
ans<-get_raw_cop(nametag=nametag,dataloc=dataloc,resloc=resloc)
saveRDS(ans,paste(resloc,nametag,"_coplist.RDS",sep=""))

dataloc<-"./Data/BMR/simulations/sims/SomewhatLeftTailDep/"
nametag<-"SomewhatL"
ans<-get_raw_cop(nametag=nametag,dataloc=dataloc,resloc=resloc)
saveRDS(ans,paste(resloc,nametag,"_coplist.RDS",sep=""))

dataloc<-"./Data/BMR/simulations/sims/SomewhatRightTailDep/"
nametag<-"SomewhatR"
ans<-get_raw_cop(nametag=nametag,dataloc=dataloc,resloc=resloc)
saveRDS(ans,paste(resloc,nametag,"_coplist.RDS",sep=""))

#========================================= Now do a NPA for each data==========================================

do_evo_NPA<-function(coplist,bound){
  
  # initiate to store the results from 100 copulas
  Corl_list<-c()
  Coru_list<-c()
  Pl_list<-c()
  Pu_list<-c()
  D2u_list<-c()
  D2l_list<-c()
  CorlmCoru_list<-c()
  PlmPu_list<-c()
  D2umD2l_list<-c()
  
  for(i in 1:100){
    
    v<-coplist[[i]]
    
    #Cor stat
    Corl<-Corbds(vi=v[,1],vj=v[,2],lb=0,ub=bound)
    Coru<-Corbds(vi=v[,1],vj=v[,2],lb=(1-bound),ub=1)
    CorlmCoru<-Corl-Coru
    
    #P stat
    res_l<-Pbds(vi=v[,1],vj=v[,2],lb=0,ub=bound)
    Pl<-res_l$abs_res
    
    res_u<-Pbds(vi=v[,1],vj=v[,2],lb=(1-bound),ub=1)
    Pu<-res_u$abs_res
    PlmPu<-Pl-Pu
    
    
    #D2 stat
    D2l<-D2bds(vi=v[,1],vj=v[,2],lb=0,ub=bound)
    D2u<-D2bds(vi=v[,1],vj=v[,2],lb=(1-bound),ub=1)
    D2umD2l<-D2u-D2l
    
    Corl_list<-c(Corl_list,Corl)
    Coru_list<-c(Coru_list,Coru)
    Pl_list<-c(Pl_list,Pl)
    Pu_list<-c(Pu_list,Pu)
    D2u_list<-c(D2u_list,D2u)
    D2l_list<-c(D2l_list,D2l)
    CorlmCoru_list<-c(CorlmCoru_list,CorlmCoru)
    PlmPu_list<-c(PlmPu_list,PlmPu)
    D2umD2l_list<-c(D2umD2l_list,D2umD2l)
    
  }
  
  res<-list(Corl_list=Corl_list,
            Coru_list=Coru_list,
            Pl_list=Pl_list,
            Pu_list=Pu_list,
            D2u_list=D2u_list,
            D2l_list=D2l_list,
            CorlmCoru_list=CorlmCoru_list,
            PlmPu_list=PlmPu_list,
            D2umD2l_list=D2umD2l_list
            )
    
  return(res)
  
}

#---------------Now call the function---------------------------------

nametag<-"BivNormal"
coplist<-readRDS("./Results/BMR_results/evo_cop/BivNormal_coplist.RDS")
resloc<-"./Results/BMR_results/evo_cop/"
ans<-do_evo_NPA(coplist=coplist,bound = 0.2)
saveRDS(ans,paste(resloc,nametag,"_NPA_results.RDS",sep=""))

pdf(paste(resloc,nametag,"_hist_taildep.pdf",sep=""),height=2,width=6)
op<-par(mfrow=c(1,3),mar=c(5,5,0.5,1),mgp=c(3.5,1,0))

res<-ans$CorlmCoru_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(cor[l]-cor[u]),cex.lab=2,cex.axis=2,main="")
     #main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$PlmPu_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(P[l]-P[u]),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$D2umD2l_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(D[u]^2 - D[l]^2),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

par(op)
dev.off()

#----------------------
nametag<-"ExtremeL"
coplist<-readRDS("./Results/BMR_results/evo_cop/ExtremeL_coplist.RDS")
ans<-do_evo_NPA(coplist=coplist,bound = 0.2)
resloc<-"./Results/BMR_results/evo_cop/"
saveRDS(ans,paste(resloc,nametag,"_NPA_results.RDS",sep=""))

pdf(paste(resloc,nametag,"_hist_taildep.pdf",sep=""),height=2,width=6)
op<-par(mfrow=c(1,3),mar=c(5,5,0.5,1),mgp=c(3.5,1,0))

res<-ans$CorlmCoru_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(cor[l]-cor[u]),cex.lab=2,cex.axis=2,main="")
#main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$PlmPu_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(P[l]-P[u]),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$D2umD2l_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(D[u]^2 - D[l]^2),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

par(op)
dev.off()

#-------------------------
nametag<-"ExtremeR"
coplist<-readRDS("./Results/BMR_results/evo_cop/ExtremeR_coplist.RDS")
ans<-do_evo_NPA(coplist=coplist,bound = 0.2)
resloc<-"./Results/BMR_results/evo_cop/"
saveRDS(ans,paste(resloc,nametag,"_NPA_results.RDS",sep=""))

pdf(paste(resloc,nametag,"_hist_taildep.pdf",sep=""),height=2,width=6)
op<-par(mfrow=c(1,3),mar=c(5,5,0.5,1),mgp=c(3.5,1,0))

res<-ans$CorlmCoru_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(cor[l]-cor[u]),cex.lab=2,cex.axis=2,main="")
#main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$PlmPu_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(P[l]-P[u]),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$D2umD2l_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(D[u]^2 - D[l]^2),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

par(op)
dev.off()

#------------------------
nametag<-"SomewhatL"
coplist<-readRDS("./Results/BMR_results/evo_cop/SomewhatL_coplist.RDS")
ans<-do_evo_NPA(coplist=coplist,bound = 0.2)
resloc<-"./Results/BMR_results/evo_cop/"
saveRDS(ans,paste(resloc,nametag,"_NPA_results.RDS",sep=""))

pdf(paste(resloc,nametag,"_hist_taildep.pdf",sep=""),height=2,width=6)
op<-par(mfrow=c(1,3),mar=c(5,5,0.5,1),mgp=c(3.5,1,0))

res<-ans$CorlmCoru_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(cor[l]-cor[u]),cex.lab=2,cex.axis=2,main="")
#main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$PlmPu_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(P[l]-P[u]),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$D2umD2l_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(D[u]^2 - D[l]^2),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

par(op)
dev.off()

#----------------------------
nametag<-"SomewhatR"
coplist<-readRDS("./Results/BMR_results/evo_cop/SomewhatR_coplist.RDS")
ans<-do_evo_NPA(coplist=coplist,bound = 0.2)
resloc<-"./Results/BMR_results/evo_cop/"
saveRDS(ans,paste(resloc,nametag,"_NPA_results.RDS",sep=""))

pdf(paste(resloc,nametag,"_hist_taildep.pdf",sep=""),height=2,width=6)
op<-par(mfrow=c(1,3),mar=c(5,5,0.5,1),mgp=c(3.5,1,0))

res<-ans$CorlmCoru_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(cor[l]-cor[u]),cex.lab=2,cex.axis=2,main="")
#main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$PlmPu_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(P[l]-P[u]),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

res<-ans$D2umD2l_list
mres<-round(median(res),4)
hist(res,border=F,col="grey",breaks = 10,xlab=expression(D[u]^2 - D[l]^2),cex.lab=2,cex.axis=2,main="")
#     main=paste("median=",mres,sep=""))
abline(v=0,lty="dashed")

par(op)
dev.off()














