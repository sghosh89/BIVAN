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
#results of D2bds(vi,vj,0,.5) and D2bds(vi,vj,0.5,1).
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

#Calculates the portion of the Spearman correlation 
#that is due to points in from a copula that are
#between two bounds.
#
#Args
#vi, vj       Coordinates of points from a copula
#lb, ub       Lower and upper bounds between 0 and 1
#
#Output
#The portion of the Spearman correlation.
#
#Examples
#The earlier function CorlCoru should give the combined
#results of Corbds(vi,vj,0,.5) and Corbds(vi,vj,0.5,1).
#
Corbds<-function(vi,vj,lb,ub)
{
  #get mean and variance
  vi_mean<-mean(vi)
  vj_mean<-mean(vj)
  var_vi<-var(vi)
  var_vj<-var(vj)
  
  #compute the indices of the points between the bounds
  inds<-which(vi+vj>2*lb & vi+vj<2*ub)
  
  #get the portion of the Spearman
  res<-sum((vi[inds]-vi_mean)*(vj[inds]-vj_mean))/((length(vi)-1)*sqrt(var_vi*var_vj))

  return(res)  
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