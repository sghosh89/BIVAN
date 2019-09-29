# This function gives you a list of two matrices where its [i,j] element are the AICw and AIC 
# for a fitted normal copula for the location pair [i,j]
# Input : 
#       r is the RDS file which came as an output from FittingCopula_selective_loc.R

normal_cop_AICw_matrix<-function(r){
  normal_cop_AICw<-r$gfc_numBS
  normal_cop_AIC<-normal_cop_AICw
  #matrix(NA, nrow=dim(r$gfc_p_CvM)[1], ncol=dim(r$gfc_p_CvM)[2])
  #colnames(normal_cop_AICw)<-colnames(r$gfc_p_CvM)
  #rownames(normal_cop_AICw)<-rownames(r$gfc_p_CvM)
  
  for(i in c(1:dim(normal_cop_AICw)[1])){
    for(j in c(1:dim(normal_cop_AICw)[1])){
      if(is.finite(r$gfc_numBS[i,j])){
        x<-r$info_ord_copcode[i,j,]
        norm_ind<-which(x==1)
        norm_AICw<-r$info_ord_AICw[i,j,norm_ind]
        normal_cop_AICw[i,j]<-norm_AICw
        norm_AIC<-r$info_ord_AIC[i,j,norm_ind]
        normal_cop_AIC[i,j]<-norm_AIC
      }
    }
  } 
  
  return(list(normal_cop_AICw=normal_cop_AICw,
              normal_cop_AIC=normal_cop_AIC))
  
}

#--------------
#check the function
#r<-readRDS("Plankton_North_Sea_CopulaFit_selecloc_species_16.RDS")
#z<-normal_cop_AICw_matrix(r=r)








