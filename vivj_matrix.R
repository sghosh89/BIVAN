# Select any two location pair [i,j] for a copula of species sp
# This function gives you a matrix with vi and vj as two columns
library(VineCopula)
vivj_matrix<-function(d_allsp,sp,i,j){
  
  ds1<-d_allsp[[sp]][[i]]
  ds2<-d_allsp[[sp]][[j]]
  #----------------------------
  a1<-ds1$Year[1]
  a2<-ds2$Year[1]
  a3<-ds1$Year[dim(ds1)[1]]
  a4<-ds2$Year[dim(ds2)[1]]
  year_s<-max(a1,a2)
  year_e<-min(a3,a4)
  ind_s1<-which(ds1$Year==year_s)
  ind_s2<-which(ds2$Year==year_s)
  ind_e1<-which(ds1$Year==year_e)
  ind_e2<-which(ds2$Year==year_e)
  ds1<-ds1[ind_s1:ind_e1,]
  ds2<-ds2[ind_s2:ind_e2,]
  # Omitting the years and data containing NA in either d1 or d2 
  #from both d1 and d2
  if(anyNA(ds1$Dat)==T | anyNA(ds2$Dat)==T){
    ind_na1<-which(is.na(ds1$Dat))
    ind_na2<-c(ind_na1,which(is.na(ds2$Dat)))
    ind_na<-unique(ind_na2)
    
    d1Dat<-ds1$Dat[-ind_na]
    d2Dat<-ds2$Dat[-ind_na]
    Years<-ds1$Year[-ind_na]
    d1<-data.frame(Year=Years,Dat=d1Dat)
    d2<-data.frame(Year=Years,Dat=d2Dat)
  } else {
    d1<-ds1
    d2<-ds2
  }
  #get ranks modified now
  vi<-VineCopula::pobs(d1$Dat)
  vj<-VineCopula::pobs(d2$Dat)
  Years<-d1$Year
  #-------------------------
  #n_datapt<-length(vi)
  #--------------------
  #plot(vi,vj,type="p")
  #-------------------------
  mat<-as.matrix(cbind(vi,vj))
  return(mat)  
}


