require(copula)

#This function takes two time series measured at the same times 
#in different locations and creates surrogate datasets which
#are statistically similar except the copula has been changed 
#to any specified copula. Any one-parameter 2d copula implemented 
#in the copula package can be used.
#
#Args
#m          A T by 2 matrix, where T is the length of the time 
#             series. Assumed no NAs.
#targetcop  The target copula. This is a copula object from the 
#             copula package.
#numsurrog  The desired number of surrogate datasets
#
#Output
#A T by 2 by numsurrog array, surrogs. The ith surrogate is stored
#surrog[,,i].
#
copsurrog2d<-function(m,targetcop,numsurrog)
{
  T<-dim(m)[1]
  
  #get the kendall correlation of the data
  kcor<-cor(m[,1],m[,2],use="pairwise.complete.obs",method="kendall")
  
  #Get the parameter of targetcop that produces this same
  #kendall correlation and change the parameter to that.
  #This won't work if targetcop is a copula family with more 
  #than 1 parameter.
  targetcop@parameters<-iTau(targetcop,kcor)
  
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

#***basic tests of the function
#datcop<-claytonCopula(5,2)
#numpts<-1000
#m<-rCopula(numpts,datcop)
#plot(m[,1],m[,2],type='p')

#tarcop<-gumbelCopula(3,2)

#res<-copsurrog2d(m,tarcop,3)
#dim(res)

#plot(res[,1,1],res[,2,1],type='p')
#BiCopGofTest(res[,1,1],res[,2,1],family=4)
#BiCopGofTest(res[,1,2],res[,2,2],family=4)
#BiCopGofTest(res[,1,3],res[,2,3],family=4)

#sum(sort(res[,1,1])==sort(m[,1]))
#sum(sort(res[,2,1])==sort(m[,2]))
#sum(sort(res[,1,2])==sort(m[,1]))
#sum(sort(res[,2,2])==sort(m[,2]))
#sum(sort(res[,1,3])==sort(m[,1]))
#sum(sort(res[,2,3])==sort(m[,2]))

#***test to see whether changing the copula
#markedly influences the distribution
#of the mean

#make the data
#datcop<-rotCopula(claytonCopula(5,2))
#d<-rCopula(10000,datcop)
#d<-qnorm(d,mean=0,sd=1)

#change the copula
#tarcop<-claytonCopula(3,2)
#sur<-copsurrog2d(d,tarcop,1000)

#see how the skewness of the spatial average 
#is affected
#source("SkewnessAnd3CentMom.R")
#skdat<-myskns(apply(FUN=mean,X=d,MARGIN=1))
#mnres<-apply(FUN=mean,X=sur,MARGIN=c(1,3))
#sksur<-apply(FUN=myskns,X=mnres,MARGIN=2)
#sum(sksur<skdat)/length(sksur)
#range(sksur) #all negative, consistent with understanding
#skdat #positive, also consistent
