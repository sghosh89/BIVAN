source("getcopula.R")
#--------------------------------------------------------------------------------------
# This function (mainly written for phylodata for contrasts : containing path length and different traits)
# takes a dataframe and remove the 1st column of that df and then calculates stuffs on that 
# modified dataframe to return output

# Note : to get copula between traits use usual data cleaning but to find copula between pathlength and trait 
# take absolute value of trait only

# Input :
#        data : a dataframe 
#        i, j : column indices of data frame to get bivariate data set
#        ploton : logical for option plots of raw data and it's copula representation
# Output : 
#        a cop rank matrix for the raw data
phylo_dat_fn<-function(data,i,j,ploton=T){
  
  data<-data[,2:dim(data)[2]]
  
  x<-c(i,j)
  if(1 %in% x){
    ind<-which(x!=1)
    data[,x[ind]]<-abs(data[,x[ind]]) # replace by its absolute value if it's comparison btw pathlength & traits
    d<-data[,c(i,j)] 
  }else{
    d<-data[,c(i,j)] 
  }

  v_phylo<-getcopula(d=d,rankon = T,ploton = ploton)
  
  #if(ploton==T){
    #op<-par(mfrow=c(1,2))
    #plot(d[,1],d[,2],col="blue")
    #getcopula(d=d,rankon = T,ploton = T)
    #par(op)
  #}
  return(v_phylo)
}




