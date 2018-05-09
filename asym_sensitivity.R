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
# Check the function
#N<-1000
#rho<-0.8
#del<-1 #???
#alpha<-0.5
  
#sigma_glob<-matrix(c(1,rho,rho,1),2,2)
#sigma_loc<-matrix(c(del,0,0,del),2,2)

#mg<-rmvnorm(N,mean=c(0,0),sigma=sigma_glob)
#ml<-rmvnorm(N,mean=c(0,0),sigma=sigma_loc)
  
#g_noise_inc<-apply(mg,MARGIN = 2,FUN = glob_noise_fn, decrease=F)
#g_noise_dec<-apply(mg,MARGIN = 2,FUN = glob_noise_fn, decrease=T)
#h_noise<-alpha*ml
#-----------------------------------------------------------------------
# This function generates global and local noise matrices 
#            which will be called as input in next simulator function
# Input:      
#       N : number of global and local noise points you want
#       rho : covariance of global noise matrix
# Output:
#       a list of two matrices (mg for global noise mat) and (ml for local noise matrix)
get_noise_glob_loc_mat<-function(N,rho){
  
  sigma_glob<-matrix(c(1,rho,rho,1),2,2)
  sigma_loc<-matrix(c(1,0,0,1),2,2)
  
  mg<-rmvnorm(N,mean=c(0,0),sigma=sigma_glob)
  ml<-rmvnorm(N,mean=c(0,0),sigma=sigma_loc)
  
  return(list(mg=mg,
              ml=ml))
  
}
#----------------------------------------------------------------------------------
# This function simulates the equation to test asym sensitivity :
#      P(i,t+1)=beta*P(i,t)+globalnoise(i,t)+alpha*localnoise(i,t)
# Input :
#       mg : global noise matrix (dim = N by 2)
#       ml : local noise matrix (dim = N by 2)
#       rho : covariance of global noise matrix
#       alpha : a number : local noise coeeficient
#       beta : growth factor of pop
#       decrease : logical (controls the nature of global noise)
#       numkeep_last : number of points you want to keep for the pop matrix from the end
# Output :
#     pop matrix with dim. = numkeep_last by 2

Simulator_asym_sens<-function(mg,ml,alpha,beta,p0=c(0,0),decrease,numkeep_last){
  
  N<-nrow(mg)
  g_noise<-apply(mg,MARGIN = 2,FUN = glob_noise_fn, decrease=decrease)
  h_noise<-alpha*ml
  
  res<-matrix(NA,N+1,2)
  res[1,]<-p0
  
  for (counter in 2:(N+1)){
    res[counter,]<-(res[counter-1,]*beta)+g_noise[counter-1,]+h_noise[counter-1,]
  }
  
  pop<-tail(res,numkeep_last)
  rownames(pop)<-c()
  
  pop<-pobs(pop)
  
  return(pop)
  
} 
#--------------------------------------------------------------------------------------




