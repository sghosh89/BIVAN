require(stats)

#This function compares the empirical cdf of a univariate dataset to a cloud
#of empirical cdfs of surrogate datasets of the same format, visually.
#
#Args
#d          A numeric vector, assumed no NAs
#m          A numeric matrix, number of rows equal to the length of d. Assumed
#             no NAs.
#
#Output: a list with these named elements
#dcdf       An ecdf object (see the ecdf function in the stats package)
#x          A vector of x axis values spanning from min(d) to max(d)
#mcdf       The cdfs of the columns of m evaluated at the values of x
#prob       Quantiles to use in plotting
#
#A plot is generated showing the quantiles of the cloud of empirical cdfs of 
#the columns of m, with the empirical cdf of d plotted on top in red.
#
#A second plot is generated showing the fraction of surrogate empirical cdfs
#the real empirical cdf is greater than.
#
empcdfcloud<-function(d,m,probs=c(.005,.5,.995))
{
  #get the empirical cdf of d and plot
  dcdf<-ecdf(d)
  plot(dcdf,main='',col='red')
  
  #compute the surrogate empirical cdfs
  x<-seq(from=min(d),to=max(d),length.out=250)
  mcdf<-matrix(NA,length(x),dim(m)[2])
  for (counter in 1:length(x))
  {
    mcdf[counter,]<-apply(FUN=sum,MARGIN=2,X=(m<=x[counter]))/length(d)
  }
  qcdf<-apply(FUN=quantile,X=mcdf,MARGIN=1,probs=probs)
  
  #add the quantiles to the plot
  for (counter in 1:(dim(qcdf)[1]))
  {
    lines(x,qcdf[counter,],type='l',lty='solid',col='grey')
  }

  #make the second plot
  dcdfy<-dcdf(x)
  py<-apply(FUN=sum,X=(matrix(dcdfy,length(dcdfy),dim(mcdf)[2])>mcdf),MARGIN=1)/dim(mcdf)[2]
  plot(x,py,type='l',xlab='x',ylab='Frac surrog cdfs less than data cdf',ylim=c(0,1))
  lines(range(x),rep(.025,2),lty='dashed')
  lines(range(x),rep(.975,2),lty='dashed')
  
  return(list(dcdf=dcdf,x=x,mcdf=mcdf))
}

##test
#d<-rnorm(1000)
#m<-matrix(rnorm(1000*100,mean=.1),1000,100)
#res<-empcdfcloud(d,m)
#
#m<-matrix(rnorm(1000*100,mean=.5),1000,100)
#res<-empcdfcloud(d,m)
#
#m<-matrix(rnorm(1000*100,mean=0,sd=1.3),1000,100)
#res<-empcdfcloud(d,m)
