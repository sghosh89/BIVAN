#-------------------------
library(VineCopula)
library(mvtnorm)
#----------------------------
# This function takes a vector (x) and a logical tag and then transform them accordingly
glob_noise_fn<-function(x,decrease){
  
  if(decrease==T){
    x<--x
  }
  
  for(i in 1:length(x)){
    if(x[i]<0){
      x[i]<-x[i]
    }else{
      x[i]<-0
    }
  }
  
  if(decrease==T){
    stopifnot(all(-x>=0))  # This is a check
    return(-x)
  }else{
    stopifnot(all(x<=0))   # This is a check
    return(x)
  }
}
#-----------------------------------------------------------------------
# This function generates global and local noise matrices 
#            which will be called as input in next simulator function
# Input:      
#       N : number of global and local noise points you want
#       r : covariance of global noise matrix
# Output:
#       a list of two matrices (mg for global noise mat) and (ml for local noise matrix)
get_noise_glob_loc_mat<-function(N,r){
  
  sigma_glob<-matrix(c(1,r,r,1),2,2)
  sigma_loc<-matrix(c(1,0,0,1),2,2)
  
  mg<-rmvnorm(N,mean=c(0,0),sigma=sigma_glob)
  ml<-rmvnorm(N,mean=c(0,0),sigma=sigma_loc)
  
  return(list(mg=mg,
              ml=ml))
  
}
#----------------------------------------------------------------------------------
# This function simulates the equation to test asym sensitivity :
#      P(i,t+1)=b*P(i,t)+globalnoise(i,t)+a*localnoise(i,t)
# Input :
#       mg : global noise matrix (dim = N by 2)
#       ml : local noise matrix (dim = N by 2)
#       r : covariance of global noise matrix
#       a : a number : local noise coeeficient
#       b : growth factor of pop
#       decrease : logical (controls the nature of global noise)
#       numkeep_last : number of points you want to keep for the pop matrix from the end
# Output :
#     A list of two :
#     pop matrix with dim. = numkeep_last by 2
#     g_noise : global noise

Simulator_asym_sens<-function(mg,ml,a,b,p0=c(0,0),decrease,numkeep_last){
  
  N<-nrow(mg)
  g_noise<-apply(mg,MARGIN = 2,FUN = glob_noise_fn, decrease=decrease)
  h_noise<-a*ml
  
  res<-matrix(NA,N+1,2)
  res[1,]<-p0
  
  for (counter in 2:(N+1)){
    res[counter,]<-(res[counter-1,]*b)+g_noise[counter-1,]+h_noise[counter-1,]
  }
  
  pop<-tail(res,numkeep_last)
  rownames(pop)<-c()
  
  pop<-pobs(pop)
  
  return(list(pop=pop,
              g_noise=g_noise))
  
} 
#--------------------------------------------------------------------------------------

