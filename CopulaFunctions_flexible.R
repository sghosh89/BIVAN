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
Corbds<-function(vi,vj,lb,lu)
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
#Outputs
#
Pbds<-function(vi,vj,lb,ub){
  
  if(lb>=ub){cat("error : lb is greater or equal to ub","\n")}
  
  # when two boundary lines are on the right side of vi+vj=1 line  
  if((2*lb>=1) & (2*ub >=1)){
    d_max<-lb*sqrt(2)
    a<-2*sqrt(2)*(ub-1)
    b<-2*d_max
    h<-abs(2*(ub-lb))/sqrt(2)
    deno<-2*(((lb-1)^2)-((ub-1)^2))
    
    dist_Si1<-c()
    Si1<-c()
    for (di in seq(from=0, to=(0.5*a), by=(0.5*a)/1000)){ 
      dist_Si1<-c(dist_Si1,di)
      a1<-2*di*h
      Si1<-c(Si1,a1) 
    } 
    
    dist_Si2<-c()
    Si2<-c()
    ditemp<-d_max-(a/2)
    for(di in seq(from=0, to=ditemp, by=ditemp/1000)){
      dist_Si2<-c(dist_Si2,(tail(dist_Si1,1)+di))
      a2<-2*(0.5*di*(h+(h-di)))
      a2<-tail(Si1,1)+a2
      Si2<-c(Si2,a2)
    }  
    
    dist_Si<-c(dist_Si1,dist_Si2[-1])
    Si<-c(Si1,Si2[-1])/deno
    
  # when two boundary lines are on the left side of vi+vj=1 line  
  }else if((2*lb<=1) & (2*ub <=1)){
    d_max<-ub*sqrt(2)
    a<-2*sqrt(2)*lb
    b<-2*d_max
    h<-abs(2*(ub-lb))/sqrt(2)
    deno<-2*((ub^2)-(lb^2))
    
    dist_Si1<-c()
    Si1<-c()
    for (di in seq(from=0, to=(0.5*a), by=(0.5*a)/1000)){ 
      dist_Si1<-c(dist_Si1,di)
      a1<-2*di*h
      Si1<-c(Si1,a1) 
    } 
    
    dist_Si2<-c()
    Si2<-c()
    ditemp<-d_max-(a/2)
    for(di in seq(from=0, to=ditemp, by=ditemp/1000)){
      dist_Si2<-c(dist_Si2,(tail(dist_Si1,1)+di))
      a2<-2*(0.5*di*(h+(h-di)))
      a2<-tail(Si1,1)+a2
      Si2<-c(Si2,a2)
    }  
    
    dist_Si<-c(dist_Si1,dist_Si2[-1])
    Si<-c(Si1,Si2[-1])/deno
    
  # when lower boundary lines are on the left side of vi+vj=1 line and upper on the other side    
  }else{
    d_max<-1/sqrt(2)
    lblen<-sqrt(2)*lb
    ublen<-abs(sqrt(2)*(ub-1))
    di1<-min(lblen,ublen)
    di2<-max(lblen,ublen)
    h<-abs(2*(ub-lb))/sqrt(2)
    deno<-(0.5-(2*(lb^2)))+(0.5-(2*((ub-1)^2)))
    
    dist_Si1<-c()
    Si1<-c()
    for (di in seq(from=0, to=di1, by=di1/1000)){ 
      dist_Si1<-c(dist_Si1,di)
      a1<-2*di*h
      Si1<-c(Si1,a1) 
    } 
    
    dist_Si2<-c()
    Si2<-c()
    ditemp1<-di2-di1
    for(di in seq(from=0, to=ditemp1, by=ditemp1/1000)){
      dist_Si2<-c(dist_Si2,(tail(dist_Si1,1)+di))
      a2<-2*(0.5*di*(h+(h-di)))
      a2<-tail(Si1,1)+a2
      Si2<-c(Si2,a2)
    }
    
    base<-(h-ditemp1)
    
    dist_Si3<-c()
    Si3<-c()
    ditemp2<-d_max-di2
    for(di in seq(from=0, to=ditemp2, by=ditemp2/1000)){
      dist_Si3<-c(dist_Si3,(tail(dist_Si2,1)+di))
      a3<-2*(0.5*di*(base+(base-(2*di))))
      a3<-tail(Si2,1)+a3
      Si3<-c(Si3,a3)
    }  
    
    dist_Si<-c(dist_Si1,dist_Si2[-1],dist_Si3[-1])
    Si<-c(Si1,Si2[-1],Si3[-1])/deno
  
  }
  
  #---------------------------------------------------------
  #Au_Si<-0
  #for(i in 1:(length(Si)-1)){
  #  Au_Si<-Au_Si+(Si[i]*(dist_Si[i+1]-dist_Si[i]))
  #}
  
  # analytical ingration
  
  # when lower boundary lines are on the left side of vi+vj=1 line and upper on the other side 
  if((2*lb<1) & (2*ub >1)){
    Au_Si<-(ditemp1^2)*(h-base)
    Au_Si<-Au_Si+((di1^3)/3)
    Au_Si<-Au_Si+(base*(ditemp2^2))
    Au_Si<-Au_Si-((ditemp2^3)/3)
    Au_Si<-Au_Si/deno
  
  }else{
    # when both boundaries are on the same side of vi+vj=1 line
    Au_Si<-0.25*h*(a^2)+(h*(d_max^2))+((a^3)/24)-((d_max^3)/3)  
    Au_Si<-Au_Si/deno
    
  }
  
  #---------------------------------------------------------
  
  inds<-which(vi+vj>2*lb & vi+vj<2*ub)
  dist_sort<-sort(abs(vi[inds]-vj[inds])/sqrt(2))
  
  dpt_c<-0
  dpt_uniq<-c(0)
  dist_sort_df<-as.data.frame(table(dist_sort))
  for (i in 1:length(dist_sort_df$Freq)){
    dpt_c<-dpt_c+(dist_sort_df$Freq[i]/length(dist_sort))
    dpt_uniq<-c(dpt_uniq,dpt_c)
  }
  dpt_uniq<-c(dpt_uniq,1)
  dist_dpt_uniq<-c(0,as.numeric(as.character(dist_sort_df$dist_sort)),d_max)
  
  nrep<-2
  S<-rep(dpt_uniq,each=nrep)
  dist_S<-rep(dist_dpt_uniq[2:length(dist_dpt_uniq)],each=nrep)
  dist_S<-c(0,dist_S,d_max)
  #plot(dist_S,S,type="l",col="red")
  
  #integrals under the functions S
  Au_S<-0
  
  for(i in 1:(length(dpt_uniq)-1)){
    Au_S<-Au_S+(dpt_uniq[i]*(dist_dpt_uniq[i+1]-dist_dpt_uniq[i]))
  }
  #-----------------------------------------------------------------------
  res<-Au_S-Au_Si
  
  return(list(dist_S=dist_S,
              S=S,
              dist_Si=dist_Si,
              Si=Si,
              abs_res=abs(res)))
 
}



