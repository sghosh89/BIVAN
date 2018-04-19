library(VineCopula)

#This function takes a bivariable dataset (matrix or dataframe) and gets the
#associated sample from the underlying copula by taking ranks.
#
#Args
#d          An N by 2 matrix or data frame
#rankon     Should ranks be taken? This could be F, for instance, if data 
#             were simulated from a copula.
#ploton     Should a plot of the copula be generated?
#
#Note: An NA means that whole row is automatically removed.
#
#Output : a plot(optional) and rankings as a matrix
#
getcopula<-function(d,rankon=T,ploton=F)
{
  #get rid of NAs
  inds<-which(is.finite(d[,1]) & is.finite(d[,2]))
  d<-d[inds,]
  
  #get ranks
  if (rankon)
  {
    v<-VineCopula::pobs(d)
  }else
  {
    v<-d
  }
  
  #plot
  if(ploton==T){
    op<-par(mar=c(5.1, 4.1, 0.2, 2.1))
    plot(v[,1],v[,2],type="p",col="red",xlab="u",ylab="v",
         cex.lab=2,cex.axis=1.5,
         xlim=c(0,1),ylim=c(0,1),asp=1)
    rect(0,0,1,1)
    lines(c(0,1),c(0,1),type='l')
    par(op)
  }
  
  return(v)
}
