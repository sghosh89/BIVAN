# Function to generate a table for npa results-summary with bivariate data, 
# this one about asymmetry of tail dependence

fractxt<-function(fracval,numsurrog){
  ptxt<-''
  if (fracval>.5*numsurrog){                        
    ptxt<-paste(">",fracval,sep='')
  } 
  if (fracval<=.5*numsurrog){
    ptxt<-paste("<",numsurrog-fracval,sep='')
  }
  return(ptxt)
}

make_table_stat_fn<-function(data){
  
  #atdres<-data.frame(Statistic=c("$\\cor_{0,0.1}-\\cor_{0.9,1}$",
  #                               "$\\Ps_{0,0.1}-\\Ps_{0.9,1}$",
  #                               "$\\Dsq_{0.9,1}-\\Dsq_{0,0.1}$"),
  atdres<-data.frame(Statistic=c("cor$_{FS}$-cor$_{LS}$",
                                 "P$_{FS}$-P$_{LS}$",
                                 "D$_{LS}^2$-D$_{FS}^2$"),
                     Kendall=rep("",3),
                     Spearman=rep("",3),stringsAsFactors = F)
  
  atdres[1,2]<-fractxt(data$corlmcoru_frac_K,numsurrog = data$numsurrog_success_K)
  atdres[1,3]<-fractxt(data$corlmcoru_frac_S,numsurrog = data$numsurrog_success_S)
  atdres[2,2]<-fractxt(data$PlmPu_frac_K,numsurrog = data$numsurrog_success_K)
  atdres[2,3]<-fractxt(data$PlmPu_frac_S,numsurrog = data$numsurrog_success_S)
  atdres[3,2]<-fractxt(data$D2umD2l_frac_K,numsurrog = data$numsurrog_success_K)
  atdres[3,3]<-fractxt(data$D2umD2l_frac_S,numsurrog = data$numsurrog_success_S)
  
  return(atdres)
  
}
