require(copula)

#This function takes two time series measured at the same times 
#in different locations and creates surrogate datasets which
#are statistically similar except the copula has been changed 
#to any specified copula. Any one-parameter 2d copula implemented 
#in the copula package can be used as the target copula.
#
#Args
#m          A T by 2 matrix, where T is the length of the time 
#             series. Assumed no NAs.
#targetcop  The target copula. This is a copula object from the 
#             copula package.
#corpres    Should be "spearman" or "kendall" according to which
#             correlation coefficient you want to preserve.
#numsurrog  The desired number of surrogate datasets
#
#Output
#A T by 2 by numsurrog array, surrogs. The ith surrogate is stored
#surrog[,,i].
#
copsurrog2d<-function(m,targetcop,corpres,numsurrog)
{
  T<-dim(m)[1]
  
  if (corpres=="kendall")
  {
    #get the kendall correlation of the data
    kcor<-cor(m[,1],m[,2],use="pairwise.complete.obs",method="kendall")
  
    #Get the parameter of targetcop that produces this same
    #kendall correlation and change the parameter to that.
    #This won't work if targetcop is a copula family with more 
    #than 1 parameter.
    targetcop@parameters<-iTau(targetcop,kcor)
  }
  if (corpres=="spearman")
  {
    #get the spearman correlation of the data
    scor<-cor(m[,1],m[,2],use="pairwise.complete.obs",method="spearman")
    
    #Get the parameter of targetcop that produces this same
    #spearman correlation and change the parameter to that.
    #This won't work if targetcop is a copula family with more 
    #than 1 parameter.
    targetcop@parameters<-iRho(targetcop,scor)
  }
  if (corpres!="kendall" && corpres!="spearman")
  {
    stop("Error in copsurrog2d: corpres must be 'kendall' or 'spearman'")
  }
  
  #Generate a bunch of draws from targetcop in the shape of the 
  #final desired output. 
  surrogs<-array(rCopula(T*numsurrog,targetcop),
                 dim=c(T,numsurrog,2))
  surrogs<-aperm(surrogs,c(1,3,2)) #now it is T by 2 by numsurrog
  
  #The nth-smallest element of each time series surrog[,a,b] 
  #is replaced by the nth-smallest element of m[,a], for all n
  for (ca in 1:2)
  {
    om<-order(m[,ca])
    for (cb in 1:numsurrog)
    {
      os<-order(surrogs[,ca,cb])
      surrogs[os,ca,cb]<-m[om,ca]  
    }
  }
  
  return(surrogs)
}

