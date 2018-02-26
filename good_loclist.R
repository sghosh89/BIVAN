good_loclist<-function(d_allsp,sp,data_pt_thrs){
  
  lenloc<-length(d_allsp[[1]])
  
  finitedat_ind<-vector(mode = "list", length = lenloc)
  num_datapt<-c()
  for(loc in c(1:lenloc)){
    finitedat_ind[[loc]]<-which(!is.na(d_allsp[[sp]][[loc]]$Dat))
    num_datapt<-c(num_datapt,length(finitedat_ind[[loc]]))
  }
  
  good_loc<-which(num_datapt>=data_pt_thrs)
  return(good_loc)
  
}