#Functions for good estimators for 3rd central
#moment and skewness.

#Computes 3rd central moment of the data using
#the unbiased estimator. See
#http://mathworld.wolfram.com/SampleCentralMoment.html
#http://mathworld.wolfram.com/h-Statistic.html
#
#Args
#x      A numeric vector
#
#Output
#The estimated third central moment
my3cm<-function(x,na.rm=F)
{
  if (na.rm==T)
  {
    x<-x[is.finite(x)]  
  }
  
  n<-length(x)
  m3<-(1/n)*sum((x-mean(x))^3)
  h3<-n^2*m3/((n-2)*(n-1))
  return(h3)
}

#Computes skewness of the data, making use
#of the function my3cm. See
#https://en.wikipedia.org/wiki/Skewness
#
#Args
#x      A numeric vector
#
#Output
#The estimated skewness
#
myskns<-function(x,na.rm=F)
{
  if (na.rm==T)
  {
    x<-x[is.finite(x)]  
  }
  
  return(my3cm(x)/(sd(x)^3))
}
