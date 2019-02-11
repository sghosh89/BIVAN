library(copula)
library(mvtnorm)
library(matrixcalc)

#This function takes a bunch of time series measured at the same
#times in different locations and creates surrogate datasets which
#are statistically similar except the copula has been randomized
#to a copula in the normal family. "Statistcally similar" means the location
#marginals are (exactly) the same and the pairwise kendall or 
#spearman correlations of time series between locations are 
#similar.
#
#This is better than copsurrod2d because it can do multivariate data
#but it is worse that copsurrog2d because the target copula can only 
#be a normal copula.
#
#Notes: If the data have ties, Kendall cannot be used. In some cases
#the algorithm does not work because an intermediate matrix which needs
#to be positive semi-definite is not. In that case it throws an error.
#
#Args
#m          A N by n matrix, where N is the length of the time series 
#           and n is the number of time series
#corpres    Should be "spearman" or "kendall" according to which
#             correlation coefficient you want to preserve.
#numsurrog  The desired number of surrogate datasets
#
#Output
#A N by n by numsurrog array, surrogs. The ith surrogate is stored
#surrog[,,i].
#
ncsurrog<-function(m,corpres,numsurrog){
  N<-dim(m)[1]
  n<-dim(m)[2]
  
  #Some error checking
  if (corpres!="kendall" && corpres!="spearman"){
    stop("Error in ncsurrog: corpres must be 'kendall' or 'spearman'")
  }
  if (corpres=="kendall" && 
      any(apply(X=m,MARGIN=2,FUN=function(x){length(x)-length(unique(x))})>0))
  {
    stop("Error in ncsurrog: ties in time series are not allowed if corpres is kendall")
  }
  
  if (corpres=="kendall"){
    #get the kendall correlation matrix of the data
    kcor<-cor(m,use="pairwise.complete.obs",method="kendall")
  
    #get the covariance matrix of an mv normal with all
    #marginal variances 1 and with kendall correlation matrix 
    #kcor
    ncov<-iTau(normalCopula(0,2),kcor)
  }
  if (corpres=="spearman"){
    #get spearman correlation matrix of the data
    scor<-cor(m,use="pairwise.complete.obs",method="spearman")
    
    #get the covariance matrix of an mv normal with all
    #marginal variances 1 and with spearman correlation matrix 
    #scor
    ncov<-iRho(normalCopula(0,2),scor)
  }
  
  #more error checking
  if(!is.positive.semi.definite(ncov))
  {
    stop("Error in ncsurrog: ncov is not positive semidefinite")
  }
  
  #generate a bunch of mv normals in the shape of the final
  #desired output. Each row is a draw from an mv normal with 
  #mean 0 and cov matrix ncov.
  surrogs<-array(rmvnorm(N*numsurrog,mean=rep(0,n),sigma=ncov),
                 dim=c(N,numsurrog,n))
  surrogs<-aperm(surrogs,c(1,3,2)) #now it is T by n by numsurrog
  
  #The nth-smallest element of each time series surrog[,a,b] 
  #is replaced by the nth-smallest element of m[,a], for all n.
  #NAs are omitted.
  for (ca in 1:n){
    goodlocs<-is.finite(m[,ca])
    mgl<-m[goodlocs,ca]
    omgl<-order(mgl)
    for (cb in 1:numsurrog){
      sgl<-surrogs[goodlocs,ca,cb]
      osgl<-order(sgl)
      sgl[osgl]<-mgl[omgl]
      surrogs[goodlocs,ca,cb]<-sgl
      surrogs[!goodlocs,ca,cb]<-NA
    }
  }

  #***The below is what the code looked like before 
  #I added the ability handle NAs. I am keeping it around
  #just in case
  #
  #The nth-smallest element of each time series surrog[,a,b] 
  #is replaced by the nth-smallest element of m[,a], for all n
  #for (ca in 1:n)
  #{
  #  om<-order(m[,ca])
  #  for (cb in 1:numsurrog)
  #  {
  #    os<-order(surrogs[,ca,cb])
  #    surrogs[os,ca,cb]<-m[om,ca]  
  #  }
  #}
  
  return(surrogs)
}

