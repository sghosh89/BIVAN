#--------------------------------------------------------------------------
source("ExtremeTailDep.R")
#--------------------------------------------------------------------------
#The simulator of the multi-location model.
#
#Args
#p0       A length numlocs vector holding initial populations
#ns       An numsims by numlocs by numsteps array of the epsilons, 
#           where numsteps is the number of time steps you want
#
#Output
#A numsims by numlocs by numsteps+1 array of populations
#
#-------------------------------------------------------------------------------
# A function to generate dispersal matrix D (a square matrix : numlocs by numlocs)
# Input :
#       1)numlocs : number of location : an integer
#       2) d : degree of dispersal : a number in [0,1]
#       3) disp_everywhere :logical:
#             if T : gives D for linear chain model with equal dispersal everywhere
#             if F : gives D for linear chain model with equal dispersal only to nearest neighbor location
# Output : dispersal matrix D (a square matrix : numlocs by numlocs)

Dispersal_mat<-function(numlocs,d,disp_everywhere){
  D<-matrix(0,numlocs,numlocs)
  # A companion matrix (delta) that indicates how "off" a diagonal is:
  mat <- matrix(seq_len(numlocs ^ 2), ncol = numlocs)
  delta <- row(mat) - col(mat)
  
  diag(D)<-(1-d)
  
  if(disp_everywhere==T){
    D[delta!=0]<-d/(numlocs-1)
  }else{
    D[delta ==1 | delta==-1] <- d/2
    D[numlocs-1,numlocs]<-d
    D[2,1]<-d
  } 
  
  stopifnot(all(colSums(D)==1)) # this is a check
  return(D)
}
#-----------------------------------------------------------------------------------
#The simulator of the multi-location model with a dispersal matrix.
#-------Args--------------------------------------------------
#p0       A length numlocs vector holding initial populations
#ns       An numsims by numlocs by numsteps array of the epsilons, 
#           where numsteps is the number of time steps you want
#D        Dispersal matrix
#r,K,a,b  model parameters
# ext_thrs  a threshold below which populations go extinct
#model    a character specifying model name
#--------Output------------------------------------------------------
#A numsims by numlocs by numsteps+1 array of populations
#----------------------------------------------------------------------
popsim_ml_D<-function(p0,ns,D,r,K,a,b,ext_thrs,model){
  numsims<-dim(ns)[1]
  numlocs<-dim(ns)[2]
  numsteps<-dim(ns)[3]
  
  #initialization
  res<-array(NA,c(numsims,numlocs,numsteps+1))
  res[,,1]<-rep(p0,each=numsims)
  
  #step the populations forward
  for (tct in 1:numsteps){
    
    #stochastic ricker model : if r=0 then back to LC model
    
    #growth rates based on the noise
    lam_sto<-exp(ns[,,tct])  #numsims by numlocs matrix
    
    lam_ricker<-exp(r*(1-(res[,,tct]/K)))
    lam_verhulst<-1+(r*(1-res[,,tct]))
    lam_hassell<-r/((1+(a*res[,,tct]))^b)
    lam_msmith<-r/(1+((a*res[,,tct])^b))
    
    #growth prior to dispersal
    if(model=="ricker"){
      given_model<-lam_ricker
    }else if(model=="verhulst"){
      given_model<-lam_verhulst
    }else if(model=="hassell"){
      given_model<-lam_hassell
    }else if(model=="msmith"){
      given_model<-lam_msmith
    }else{
      warning("model not specified",immediate.=T,call.=T)
    }
    
    res[,,tct+1]<-given_model*lam_sto*res[,,tct] #numsims by numlocs matrix
    
    #after dispersal
    tres<-t(res[,,tct+1]) # numlocs by numsims matrix after taking transpose
    temp<-D%*%tres  # output from a matrix multiplication [(numlocs by numlocs) * (numlocs by numsims)]
    res[,,tct+1]<-t(temp) # a matrix (numsims by numlocs) 
    
    #see if there is extinction
    h<-res[,,tct+1]
    h[h<ext_thrs]<-0
    res[,,tct+1]<-h
  }
  
  return(res)
}

# function to calculate extinction risk
extrisk<-function(sims){
  totpop<-apply(FUN=sum,X=sims,MARGIN=c(1,3))
  return(apply(FUN=sum,X=(totpop==0),MARGIN=2)/dim(sims)[1])
}
#-------------------------------------------------------------------------------------------------
# function to optionally plot extinction risk agsinst time and give you back extinction risk 
#                                                            after specified numsteps 
# Input :
#       1) numsims : an integer : number of simulations
#       2) numsteps : an integer : number of time steps
#       3) numlocs : an integer : number of locations
#       4) D : dispersal matrix which is the output of Dispersal_mat function
#       5) r : model parameter : growth rate
#       6) K : Ricker model parameter : carrying capacity
#       7), 8) a, b : model parameters for Hassell model 
#       9) ext_thrs : extinction threshold below which populations go extinct
#       10) scl : a scaling factor which is multiplied with noise generated
#       11) model : a character specifying model name
#       12) ploton : logical(T or F) to get optional plot

plotter_ext_risk<-function(numsims,numsteps,numlocs,D,r,K,a,b,ext_thrs,scl,model,ploton){
 
  ns1<-retd(n=numsteps*numsims,d=numlocs,rl=1)# a righttail dep matrix(numpoints by numlocs,     
  #                                                      numpoints=numsteps*numsims)
  
  ns1<-array(ns1,c(numsteps,numsims,numlocs))# convert to an array (numsteps by numsims by numlocs)
  ns1<-aperm(ns1,c(2,3,1)) # convert to an array (numsims by numlocs by numsteps)
  
  ns1<-scl*ns1
  pops1<-popsim_ml_D(p0=rep(K,numlocs),ns=ns1,D=D,r=r,K=K,a=a,b=b,ext_thrs=ext_thrs,model=model)
  risk_right<-extrisk(pops1)
  #ps<-rep(K,numlocs)
  #print(ps)
  
  ns2<-(-ns1)
    pops2<-popsim_ml_D(p0=rep(K,numlocs),ns=ns2,D=D,r=r,K=K,a=a,b=b,ext_thrs=ext_thrs,model=model)
  risk_left<-extrisk(pops2)
  
  if(ploton==T){
    par(mfrow=c(2,1))
    plot(0:(length(risk_right)-1),risk_right,type='l',
         xlab='Time',ylab='Risk',col='blue',panel.first = grid())
    mtext("blue : right tail-dep")
    plot(0:(length(risk_left)-1),risk_left,type='l',xlab='Time',ylab='Risk',col='red',panel.first = grid())
    mtext("red : left tail-dep")
  }
  
  return(list(risk_right_after_numsteps=risk_right[numsteps+1],
              risk_left_after_numsteps=risk_left[numsteps+1]))
  
}
#----------------------------------------------------------------------------------------------------

# function to plot extinction risk agsinst degree of dispersal (d) after a certain numstep (time)

# Input :
#       1) numsims : an integer : number of simulations
#       2) numsteps : an integer : number of time steps
#       3) numlocs : an integer : number of locations
#       4) r : model parameter : growth rate
#       5) K : Ricker model parameter : carrying capacity
#       6)-7) a,b : Hassell model parameters
#       8) scl : a scaling factor which is multiplied with noise generated
#       9) ext_thrs extinction threshold below which populations go extinct
#       10) model : a character specifying model name
#       11) disp_everywhere :logical:
#             if T : gives D for linear chain model with equal dispersal everywhere
#             if F : gives D for linear chain model with equal dispersal only to nearest neighbor location
#       12) ploton : logical to get optional plot

varying_d<-function(numsims,numsteps,numlocs,r,K,a,b,ext_thrs=ext_thrs,scl,model,disp_everywhere,ploton){
  risk_right<-c()
  risk_left<-c()
  d_seq<-seq(from=0,to=1,by=0.1)
  for(d in d_seq){
    D_mat<-Dispersal_mat(numlocs=numlocs,d=d,disp_everywhere=disp_everywhere)
    riskrl<- plotter_ext_risk(numsims=numsims,numsteps = numsteps,numlocs = numlocs,D=D_mat,r=r,K=K,a=a,b=b,ext_thrs=ext_thrs,scl=scl,model=model,ploton=F)
    risk_r<-riskrl$risk_right_after_numsteps
    risk_l<-riskrl$risk_left_after_numsteps
    risk_right<-c(risk_right,risk_r)
    risk_left<-c(risk_left,risk_l)
  }
  
  if(ploton==T){
    op<-par(mfrow=c(1,2))
    plot(d_seq,risk_left,xlab='d',ylab='Risk_left',xlim=c(0,1),ylim=c(0,1),type="b",col="red",panel.first = grid())
    plot(d_seq,risk_right,xlab='d',ylab='Risk_right',xlim=c(0,1),ylim=c(0,1),type="b",col="blue",panel.first = grid())
    par(op)
    mtext(paste0("r = ", r," , numlocs = ",numlocs," , numsims = ",numsims," , numsteps = ",numsteps),side=3,line=0.2,col="navyblue")
  }
  
  return(data.frame(d_seq=d_seq,
                    risk_left=risk_left,
                    risk_right=risk_right))
  

}
#------------------------------------------------------------------------------------------------------------
