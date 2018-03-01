#Calculates the mean squared distance between copula
#points between two bounds and the main diagonal of
#the unit square.
#
#Args
#vi, vj       Coordinates of points from a copula
#lb, ub       Lower and upper bounds between 0 and 1
#
#Output
#The mean squared distance.
#
#Examples
#The earlier function D2lD2u should give the combined
#results of D2(vi,vj,0,.5) and D2(vi,vj,0.5,1).
#
D2bds<-function(vi,vj,lb,ub)
{
  inds<-which(vi+vj>2*lb & vi+vj<2*ub)

  if(length(inds)!=0)
  { 
    dsq<-0.5*(vi[inds]-vj[inds])^2
    D2<-sum(dsq)/length(dsq)
  }else
  {
    D2<-NA
  }
  
  return(D2)
}

#Dan to write an analogous Cor function

#Args
#vi, vj       Coordinates of points from a copula
#lb, ub       Lower and upper bounds between 0 and 1
#
#Output
#

Corbds<-function(vi,vj,lb,lu)
{
  
}

#Shya please write an analogous P function, fill in below.

#Args
#vi, vj       Coordinates of points from a copula
#lb, ub       Lower and upper bounds between 0 and 1
#
#Output
#
Pbds<-function(vi,vj,lb,lu)
{
  #returns a single number, computes pretty fast
}