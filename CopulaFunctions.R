#-----------------------------------------------------------------------------------------------------------------------
# THIS CODE CONTAINS ALL THE STATISTISCAL FUNCTIONS WHICH CAN BE TESTED ON THE TWO RANK VECTORS OF A GIVEN COPULA
#-----------------------------------------------------------------------------------------------------------------------

source("CopulaFunctions_flexible.R")

#--------------------------- STATISTICS 2 : correlation based Stat ---------------------------------------------------------

CorlCoru<-function(vi,vj){
  
  Corl<-Corbds(vi,vj,0,0.5)
  Coru<-Corbds(vi,vj,0.5,1)

  return(c(Corl,Coru))
}

#--------------------------- STATISTICS 4 : Stat based on the distance from right diagonal in lower and upper triangle ---------------------------------------------------------

PlPu<-function(vi,vj){
  
  ind_ll<-which(vi+vj<1)
  ind_ur<-which(vi+vj>1)
  
  if(length(ind_ll)!=0 & length(ind_ur)!=0){
    
    res_l<-Pbds(vi,vj,0,0.5)
    Pl<-res_l$abs_res
    dist_Sl_P<-res_l$dist_S
    Sl_P<-res_l$S
    dist_Si_P<-res_l$dist_Si
    Si_P<-res_l$Si
    
    res_u<-Pbds(vi,vj,0.5,1)
    Pu<-res_u$abs_res
    dist_Su_P<-res_u$dist_S
    Su_P<-res_u$S
    
    Sl_Su_Si_P<-list(dist_Sl_P=dist_Sl_P,Sl_P=Sl_P,dist_Su_P=dist_Su_P,Su_P=Su_P,dist_Si_P=dist_Si_P,Si_P=Si_P)
    
  }else{
    
    Sl_Su_Si_P<-list(dist_Sl_P=NA,Sl_P=NA,dist_Su_P=NA,Su_P=NA,dist_Si_P=NA,Si_P=NA)
    Pl<-NA
    Pu<-NA
    
  }
 
  return(list(Sl_Su_Si_P,Pl,Pu)) 
}

#--------------------------------------------------------- STATISTICS : 6 -----------------------------------------------------------------------
# get D2l : average of squared distance of points from the right diagonal of the box for lower triangle
# get D2u : average of squared distance of points from the right diagonal of the box for upper triangle

D2lD2u<-function(vi,vj)
{
  D2l<-D2bds(vi,vj,0,0.5)
  D2u<-D2bds(vi,vj,0.5,1)
  
  return(c(D2l,D2u))
}

#---------------------------------------------------------------------------------------------------------------------
#                                                     CODE ENDS HERE
#---------------------------------------------------------------------------------------------------------------------





















#and the others