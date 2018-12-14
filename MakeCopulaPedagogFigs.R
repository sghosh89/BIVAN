#' Make a pedagogical plot for a univariate copula
#'
#' @param coplist A list of 5 copula objects from the copula package
#' @param numpts Some panels show samples of this many points from the copulas
#' @param corvals A vector of five correlations for the five panels, which are
#' displayed on the plot
#' @param cortype "S", "P" or "K" for Spearman, Pearson, Kendall
#' @param thetavals A vector of five values of the copula parameter, again for display
#' @param tot_wd Total width of the figure generated, in inches
#' @param filename A filename, without extension, for saving in the working directory. A pdf 
#' is generated and saved.
#' 
#' @note For my own use at present, does little to no error catching
#' 
#' @author Daniel Reuman \email{reuman@@ku.edu}

muvcpp<-function(coplist,numpts,corvals,cortype,thetavals,tot_wd,filename)
{
  if (length(coplist)!=5)
  {
    stop("Error in muvcpp: only works with 5 copulas")
  }

  #***plot dimensions
  xmarg_ht<-.25
  ymarg_wd<-.25
  numsp<-.2
  gap<-0.1
  labgap<-.225
  titgap<-.15
  pan_wd<-(tot_wd-ymarg_wd-numsp-5*gap)/5
  if (pan_wd<=0){stop("Error in muvcpp: need a larger tot_wd")}
  pan_ht<-pan_wd
  tot_ht<-xmarg_ht+numsp+2*pan_ht+2*labgap+2*gap+titgap
  pchval<-20
  cexvalpts<-0.25
  cexvaltxt<-0.75
  cexvalpl<-1.5
  pdf(file=paste0(filename,".pdf"),width=tot_wd,height=tot_ht)
  
  x<-seq(from=0,to=1,length.out=101)
  y<-x
  dcarg<-matrix(c(rep(x,times=length(y)),rep(y,each=length(x))),length(x)*length(y),2)
  
  #top panel 1
  px<-0
  py<-1
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap+py*labgap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap+py*labgap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
  z<-dCopula(dcarg,coplist[[px+1]],log=FALSE)
  z<-matrix(z,101,101)
  contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),xaxt="n")
  axis(1,labels=FALSE)
  mtext("v",2,1.15)
  mtext("A",3,.05,cex=cexvalpl,adj=0)
  mtext(paste0(cortype,"=",round(spearvals[px+1],2),"; p=",round(thetavals[px+1],2)),3,1.25,cex=cexvaltxt)
  
  #top panel 2
  px<-1
  py<-1
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap+py*labgap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap+py*labgap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  z<-dCopula(dcarg,coplist[[px+1]],log=FALSE)
  z<-matrix(z,101,101)
  contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  mtext("B",3,.05,cex=cexvalpl,adj=0)
  mtext(paste0(cortype,"=",round(spearvals[px+1],2),"; p=",round(thetavals[px+1],2)),3,1.25,cex=cexvaltxt)
  
  #top panel 3
  px<-2
  py<-1
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap+py*labgap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap+py*labgap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
  z<-dCopula(dcarg,coplist[[px+1]],log=FALSE)
  z<-matrix(z,101,101)
  contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  mtext("C",3,.05,cex=cexvalpl,adj=0)
  mtext(paste0(cortype,"=",round(spearvals[px+1],2),"; p=",round(thetavals[px+1],2)),3,1.25,cex=cexvaltxt)
  
  #top panel 4
  px<-3
  py<-1
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap+py*labgap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap+py*labgap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)  
  z<-dCopula(dcarg,coplist[[px+1]],log=FALSE)
  z<-matrix(z,101,101)
  contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  mtext("D",3,.05,cex=cexvalpl,adj=0)
  mtext(paste0(cortype,"=",round(spearvals[px+1],2),"; p=",round(thetavals[px+1],2)),3,1.25,cex=cexvaltxt)
  
  #top panel 5
  px<-4
  py<-1
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap+py*labgap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap+py*labgap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)  
  z<-dCopula(dcarg,coplist[[px+1]],log=FALSE)
  z<-matrix(z,101,101)
  contour(x,y,log10(z),nlevels=5,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n")
  axis(1,labels=FALSE)
  axis(2,labels=FALSE)
  mtext("E",3,.05,cex=cexvalpl,adj=0)
  mtext(paste0(cortype,"=",round(spearvals[px+1],2),"; p=",round(thetavals[px+1],2)),3,1.25,cex=cexvaltxt)
  
  #bottom panel 1
  px<-0
  py<-0
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)  
  dat<-rCopula(numpts,coplist[[px+1]])
  plot(dat[,1],dat[,2],type='p',pch=pchval,cex=cexvalpts,
       xlim=c(0,1),ylim=c(0,1))
  mtext("u",1,1.15)
  mtext("v",2,1.15)
  mtext("F",3,.05,cex=cexvalpl,adj=0)
  
  #bottom panel 2
  px<-1
  py<-0
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)  
  dat<-rCopula(numpts,coplist[[px+1]])
  plot(dat[,1],dat[,2],type='p',pch=pchval,cex=cexvalpts,
       xlim=c(0,1),ylim=c(0,1),yaxt="n")
  mtext("u",1,1.15)
  axis(2,labels=FALSE)
  mtext("G",3,.05,cex=cexvalpl,adj=0)
  
  #bottom panel 3
  px<-2
  py<-0
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)    
  dat<-rCopula(numpts,coplist[[px+1]])
  plot(dat[,1],dat[,2],type='p',pch=pchval,cex=cexvalpts,
       xlim=c(0,1),ylim=c(0,1),yaxt="n")
  mtext("u",1,1.15)
  axis(2,labels=FALSE)
  mtext("H",3,.05,cex=cexvalpl,adj=0)
  
  #bottom panel 4
  px<-3
  py<-0
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)    
  dat<-rCopula(numpts,coplist[[px+1]])
  plot(dat[,1],dat[,2],type='p',pch=pchval,cex=cexvalpts,
       xlim=c(0,1),ylim=c(0,1),yaxt="n")
  mtext("u",1,1.15)
  axis(2,labels=FALSE)
  mtext("I",3,.05,cex=cexvalpl,adj=0)
  
  #bottom panel 5
  px<-4
  py<-0
  par(fig=c((ymarg_wd+numsp+px*pan_wd+px*gap)/tot_wd,
            (ymarg_wd+numsp+(px+1)*pan_wd+px*gap)/tot_wd,
            (xmarg_ht+numsp+py*pan_ht+py*gap)/tot_ht,
            (xmarg_ht+numsp+(py+1)*pan_ht+py*gap)/tot_ht),
      mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)  
  dat<-rCopula(numpts,coplist[[px+1]])
  plot(dat[,1],dat[,2],type='p',pch=pchval,cex=cexvalpts,
       xlim=c(0,1),ylim=c(0,1),yaxt="n")
  mtext("u",1,1.15)
  axis(2,labels=FALSE)
  mtext("J",3,.05,cex=cexvalpl,adj=0)
  
  dev.off()
}

#***now make some plots

numpts<-250
spearvals<-c(.1,.6,.7,.8,.9)

#normal
h<-normalCopula(.2,2,dispstr="un")
parms<-apply(X=matrix(spearvals,1,5),MARGIN=2,FUN=function(x){iRho(h,x)})
coplist<-list(normalCopula(parms[1],2,dispstr = "un"),
              normalCopula(parms[2],2,dispstr = "un"),
              normalCopula(parms[3],2,dispstr = "un"),
              normalCopula(parms[4],2,dispstr = "un"),
              normalCopula(parms[5],2,dispstr = "un"))
muvcpp(coplist,numpts,spearvals,"S",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_Normal")

#clayton
h<-claytonCopula(3,2)
parms<-apply(X=matrix(spearvals,1,5),MARGIN=2,FUN=function(x){iRho(h,x)})
coplist<-list(claytonCopula(parms[1]),
              claytonCopula(parms[2]),
              claytonCopula(parms[3]),
              claytonCopula(parms[4]),
              claytonCopula(parms[5]))
muvcpp(coplist,numpts,spearvals,"S",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_Clayton")

#survival Clayton
coplist<-lapply(FUN=rotCopula,X=coplist)
muvcpp(coplist,numpts,spearvals,"S",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_SurvClayton")

#gumbel
h<-gumbelCopula(3,2)
parms<-apply(X=matrix(spearvals,1,5),MARGIN=2,FUN=function(x){iRho(h,x)})
coplist<-list(gumbelCopula(parms[1]),
              gumbelCopula(parms[2]),
              gumbelCopula(parms[3]),
              gumbelCopula(parms[4]),
              gumbelCopula(parms[5]))
muvcpp(coplist,numpts,spearvals,"S",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_Gumbel")

#survival Gumbel
coplist<-lapply(FUN=rotCopula,X=coplist)
muvcpp(coplist,numpts,spearvals,"S",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_SurvGumbel")

#frank
h<-frankCopula(3,2)
parms<-apply(X=matrix(spearvals,1,5),MARGIN=2,FUN=function(x){iRho(h,x)})
coplist<-list(frankCopula(parms[1]),
              frankCopula(parms[2]),
              frankCopula(parms[3]),
              frankCopula(parms[4]),
              frankCopula(parms[5]))
muvcpp(coplist,numpts,spearvals,"S",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_Frank")

#joe
h<-joeCopula(3,2)
parms<-apply(X=matrix(spearvals,1,5),MARGIN=2,FUN=function(x){iTau(h,x)})
coplist<-list(joeCopula(parms[1]),
              joeCopula(parms[2]),
              joeCopula(parms[3]),
              joeCopula(parms[4]),
              joeCopula(parms[5]))
muvcpp(coplist,numpts,spearvals,"K",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_Joe")

#survival Joe
coplist<-lapply(FUN=rotCopula,X=coplist)
muvcpp(coplist,numpts,spearvals,"K",parms,tot_wd=7.5,filename="./Results/PedagogSuppMat_SurvJoe")

#' Make a pedagogical plot for a bivariate copula
#' 
#' 

#**BB1**, **BB6**, 
#**BB7**, **BB8**, **S**urvival **C**layton , **S**urvival **G**umbel, 
#**S**urvival **J**oe, **S**urvival **BB1**, **S**urvival **BB6**, **S**urvival **BB7** 
#and **S**urvival **BB8** copulae. These copulae were chosen carefully so that they covered 
#the whole range from exclusively lower (C, SG, SJ, SBB6, SBB8) or upper tail dependence 
#(SC, G, J, BB6, BB8) to both tail dependence (BB1, SBB1, BB7, SBB7) in addition to two 
#symmetric copulae (N, F). 