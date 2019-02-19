#Generates data with extreme tail dependence. Function name
#stands for "extreme tail dependence".
#
#Args
#n        Number of points to get
#d        Dimension
#rl       1 for right tail dependence, -1 for left tail dependence
#
#Output
#A matrix of dimension n by d, marginals are standard normal.
#
retd<-function(n,d,rl)
{
  m<-matrix(rnorm(n*(d+1)),n,d+1)
  flag<-(m[,1]>0) #a coin toss
  m<-m[,2:(d+1)] #an n by d matrix of iid standard normals
  
  m[flag,]<-rep(abs(m[flag,1]),times=d)
  m[!flag,]<-(-abs(m[!flag,]))
  
  if(rl==-1)
  {
    m<-(-m)  
  }
  
  return(m)
}

##test
#d<-retd(10000,2,1)
#plot(d[,1],d[,2],type='p')
#hist(d[,1],50)
#hist(d[,2],50)
#
#d<-retd(10000,2,-1)
#plot(d[,1],d[,2],type='p')
#hist(d[,1],50)
#hist(d[,2],50)
