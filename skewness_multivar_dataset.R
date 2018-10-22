#------------------------------------------------------------------------------------------------
# This function is written for measuring significant tail dependence of a multivariate dataset
#------------------------------------------------------------------------------------------------
source("ncsurrog.R")
source("SkewnessAnd3CentMom.R")
#--------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Input : ts_matrix : a multivariate dataset (a matrix with dimension = timeseries by locations)
#         loclist : a vector of selective locations 
#         numsurrog : number of replica you want to create normal cop data
#         ploton : logical
#         corpres : either "spearman" or "kendall"

# Output : a list of two : 1. a vector which measures skewness of each surrogate normal cop
#                          2. real skewness of given data
#                         and an optional histogram plot

skewness_testing<-function(ts_matrix,loclist,numsurrog,ploton,corpres){  
  
  ts_avg<-NA*numeric(length=dim(ts_matrix)[1])
  
  for(counter in c(1:dim(ts_matrix)[1])){
    temp<-mean(ts_matrix[counter,loclist],na.rm=T)
    ts_avg[counter]<-temp
  }
  
  ts_matrix_w_avg<-cbind(ts_matrix,ts_avg)
  
  skw<-apply(FUN=myskns,X=ts_matrix_w_avg,MARGIN=2)
  realskw<-tail(skw,1)
  
  #if(plot_each_hist==T){
  #  par(mfrow=c(4,3))
  # for(i in loclist){
  #   hist(ts_matrix_w_avg[,i],breaks=50,col="orange",
  #       main=paste("location",i," ,skw = ",round(skw[i],2)),xlab="ts_data",ylab="frequency")
  # }
  # hist(ts_matrix_w_avg[,dim(ts_matrix_w_avg)[2]],breaks=50,col="red",
  #      main=paste("skw = ",round(realskw,2)),xlab="spatial_avg_ts_data",ylab="frequency")
  #}
  
  # start surrogating
  #numsurrog<-5
  ts_surrog_array<-ncsurrog(m=ts_matrix[,loclist],corpres=corpres,numsurrog)
  
  arr_dim1<-dim(ts_surrog_array)[1]
  ts_surrog_avg<-matrix(NA,nrow=arr_dim1,ncol=numsurrog)
  
  for(carray in c(1:numsurrog)){
    for(jj in c(1:arr_dim1)){
      temp<-mean(ts_surrog_array[,,carray][jj,],na.rm=T)
      ts_surrog_avg[jj,carray]<-temp
    }
  }
  
  surrogskw<-apply(FUN=myskns,X=ts_surrog_avg,MARGIN=2)
  
  if(ploton==T){
    hist(surrogskw,breaks=100,main="",xlab="surrogate_skewness",ylab="frequency",cex.lab=2,cex.axis=2)
    #abline(v=realskw,col="red")
    points(x=realskw,y=0,col="red",pch=15,cex=2)
  }
  
  return(list(surrogskw=surrogskw,realskw=unname(realskw)))
  
}  

#-------------------------------------------------------------------------------------
sp_data<-function(sp,d_allsp){
  ts_matrix<-c()
  for(count in c(1:length(d_allsp[[1]]))){
    ts_matrix<-cbind(ts_matrix,d_allsp[[sp]][[count]]$Dat)
  } 
  return(ts_matrix)
}
#------------------------------------------------------------------------------------------




















