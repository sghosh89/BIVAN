library(copula)

#This function takes a bunch of time series measured at the same
#times in different locations and creates surrogate datasets. The
#location marginals of the surrogates are exactly the same as the
#input but the copula structure is that of targetcop. Kendall and
#Spearman correlations of pairs of locations are not necessarily
#the same or similar between the input m and the surrogate output
#unless targetcop has the same Kendall or Spearman found in m.
#
#This is better than copsurrog2d and ncsurrog because is can handle 
#the multivariate case, but it is worse because Kendall and Spearman
#are not preserved unless the user ensures this by choosing targetcop
#carefully. This is better than ncsurrog because any target copula
#can be used instead of just a normal copula, but worse because 
#Kendall and Spearman are not preserved (except as above).
#
#Args
#m          A N by n matrix, where N is the length of the time series 
#           and n is the number of time series
#targetcop  The target copula. This is a copula object from the 
#             copula package. Must have dimension equal to dim(m)[2].
#             The parameters as well as the family are used.
#numsurrog  The desired number of surrogate datasets
#
#Output
#A N by n by numsurrog array, surrogs. The ith surrogate is stored
#surrog[,,i].
#
copsurrognd<-function(m,targetcop,numsurrog)
{
  N<-dim(m)[1]
  n<-dim(m)[2]
  
  #error checking
  if (targetcop@dimension!=n)
  {
    stop("Error in copsurrognd: must have dim(m)[2] equal to dimension of targetcop")
  }
  
  #Generate a bunch of draws from targetcop in the shape of the 
  #final desired output. 
  surrogs<-array(rCopula(N*numsurrog,targetcop),
                 dim=c(N,numsurrog,n))
  surrogs<-aperm(surrogs,c(1,3,2)) #now it is N by n by numsurrog
  
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
  
  return(surrogs)
}