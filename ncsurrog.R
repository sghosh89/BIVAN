library(copula)
library(mvtnorm)

#This function takes a bunch of time series measured at the same
#time in different locations and creates surrogate datasets which
#are statistically similar except the copula has been randomized
#to a normal copula.
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
ncsurrog<-function(m,corpres,numsurrog)
{
  N<-dim(m)[1]
  n<-dim(m)[2]
  
  if (corpres=="kendall")
  {
    #get the kendall correlation matrix of the data
    kcor<-cor(m,use="pairwise.complete.obs",method="kendall")
  
    #get the covariance matrix of an mv normal with all
    #marginal variances 1 and with kendall correlation matrix 
    #kcor
    ncov<-iTau(normalCopula(0,2),kcor)
  }
  if (corpres=="spearman")
  {
    #get spearman correlation matrix of the data
    scor<-cor(m,use="pairwise.complete.obs",method="spearman")
    
    #get the covariance matrix of an mv normal with all
    #marginal variances 1 and with spearman correlation matrix 
    #scor
    ncov<-iRho(normalCopula(0,2),scor)
  }
  if (corpres!="kendall" && corpres!="spearman")
  {
    stop("Error in ncsurrog: corpres must be 'kendall' or 'spearman'")
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
  for (ca in 1:n)
  {
    goodlocs<-is.finite(m[,ca])
    mgl<-m[goodlocs,ca]
    omgl<-order(mgl)
    for (cb in 1:numsurrog)
    {
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

