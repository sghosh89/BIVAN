#This function gives you all the info at a selective location pair
OurBiCopSelect<-function(u1,u2,families,AICBIC="AIC",level=0.05,numBS,gofnormal){
  #first test for independence
  indeptest<-BiCopIndTest(u1,u2)$p.value
  
  # Null Hypothesis H0: Data are independent
  #if independence rejected, then get AICs for copulas and do goodness of fit stuff
  if (indeptest<level){
    #AIC stuff
    res<-data.frame(copcode=families,
                    copname=BiCopName(families, short=TRUE),
                    par1=NA,
                    par2=NA,
                    logLik=NA,
                    AIC=NA,
                    BIC=NA)
    
    for (counter in 1:(dim(res)[1])){
      tres<-BiCopEst(u1,u2,family=res[counter,1])
      res$par1[counter]<-tres$par
      res$par2[counter]<-tres$par2
      res$logLik[counter]<-tres$logLik
      res$AIC[counter]<-tres$AIC
      res$BIC[counter]<-tres$BIC
    }
    
    #gof stuff
    if (AICBIC=="AIC"){
      ind<-which.min(res$AIC)
    }
    
    if (AICBIC=="BIC"){
      ind<-which.min(res$BIC)
    }
    
    if (AICBIC!="AIC" && AICBIC!="BIC"){
      stop("Error in OurBiCopSelect: incorrect value for AICBIC")
    }
    
#    AICBIC_same_ind<-isTRUE(which.min(res$AIC)==which.min(res$BIC))
#    if(res$copcode[ind]!=2){
#      gres<-BiCopGofTest(u1,u2,family=res$copcode[ind],par=res$par1[ind],par2=res$par2[ind],method="kendall",B=numBS)
      gres<-BiCopGofTest(u1,u2,family=res$copcode[ind],method="kendall",B=numBS)
      
      gres_p_CvM<-gres$p.value.CvM
      gres_p_KS<-gres$p.value.KS
      
      gres_p_CvM_stat<-gres$statistic.CvM
      gres_p_KS_stat<-gres$statistic.KS
      
      if(gofnormal==T){
#        gres_normal_cop<-BiCopGofTest(u1,u2,family=1,par=res$par1[1],par2=res$par2[1],method="kendall",B=numBS)
        gres_normal_cop<-BiCopGofTest(u1,u2,family=1,method="kendall",B=numBS)
        gres_normal_p_CvM<-gres_normal_cop$p.value.CvM
        gres_normal_p_KS<-gres_normal_cop$p.value.KS
        
        gres_normal_p_CvM_stat<-gres_normal_cop$statistic.CvM
        gres_normal_p_KS_stat<-gres_normal_cop$statistic.KS
        
      }else{
        
        gres_normal_p_CvM<-Inf
        gres_normal_p_KS<-Inf
        
        gres_normal_p_CvM_stat<-Inf
        gres_normal_p_KS_stat<-Inf
        
      }
      
      
#    }else{
#      gres<-BiCopGofTest(u1,u2,family=res$copcode[ind],par=res$par1[ind],par2=res$par2[ind],method="white",B=numBS)
#      gres_p<-gres$p.value
      
      
#    }
    
  } else {
    
    res<-NA
    gres_p_CvM<-NA
    gres_p_KS<-NA
    
    gres_p_CvM_stat<-NA
    gres_p_KS_stat<-NA
    
    gres_normal_p_CvM<-NA
    gres_normal_p_KS<-NA
    
    gres_normal_p_CvM_stat<-NA
    gres_normal_p_KS_stat<-NA
    
  }
  
  return(list(IndepTestRes=indeptest,
              InfCritRes=res,
              GofRes_CvM=gres_p_CvM,
              GofRes_KS=gres_p_KS,
              GofRes_CvM_stat=gres_p_CvM_stat,
              GofRes_KS_stat=gres_p_KS_stat,
              GofRes_Normal_CvM=gres_normal_p_CvM,
              GofRes_Normal_KS=gres_normal_p_KS,
              GofRes_Normal_CvM_stat=gres_normal_p_CvM_stat,
              GofRes_Normal_KS_stat=gres_normal_p_KS_stat))
}











#test the above function for reasonableness
#families<-c(1,3:6,13,14,16)
#families.good<-c(1:10,13:14,16:20,104,114,204,214)
#families.mid<-c(1:10,13:14,16:20,23:24,26:30,33:34,36:40)
#families.m<-c(1:10,13:14,16:20,23:24,26:30,33:34,36:40,104,114,124,134,204,214,224,234)
#cop<-claytonCopula(5,2)
#cop<-gumbelCopula(4,2)
#cop<-frankCopula(5,2)
#dat<-rCopula(100,cop)
#allres<-OurBiCopSelect(dat[,1],dat[,2],families)
#allres

#OurBiCopSelect(u1,u2,families.good)

