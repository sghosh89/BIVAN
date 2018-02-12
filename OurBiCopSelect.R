library(VineCopula)

#This function takes a sample from a bivariate distribution
#consisting of values between 0 and 1 and:
#a) tests for independence
#b) if independence can be rejected, it fit a series of bivariate 
#   copula models and provides fitted parameters, log likelihood, 
#   AIC and BIC for each
#c) for the best copula model (by AIC or BIC), the function tests
#   for the goodness of fit in two ways
#d) optionally tests for the goodness of fit of the normal copula
#The function is a re-implementation of the BiCopSelect function
#in the VineCopula package, re-done to ensure we could control
#what we were getting.
#
#Args
#u1, u2       The sample. These can be obtained, for instance, 
#               by applying pobs in the VineCopula package to
#               samples from any bivariate distribution or to 
#               any bivariate dataset
#level        Significance level to be used for the independence
#               test
#families     The copula families to use in b above. Uses codes
#               as specified in the docs for BiCopEst in the 
#               VineCopula package.
#AICBIC       Which one to use for model selection.
#numBSsmall
#pthresh
#numBSlarge   Number of bootstraps to use for goodness of fit.
#               First, tests are run with numBSsmall, and if p-values
#               are ever less than pthresh, tests are re-run with
#               numBSlarge. Watch out, large values can mean long 
#               run times.
#gofnormal    T/F on whether a goodness of fit test for the normal
#               copula should be run.
#status       T/F - should status updates on the run be printed?
#
#Output: a list with these elements
#IndepTestRes     p-value result of the independence test
#InfCritRes       data frame of results from the information-
#                   criterion-based model fitting and selection
#GofRes_CvM       p-value result for Cramer-von-Mises-based 
#                   goodness of fit test for top-ranking copula
#                   by AIC or BIC
#GofRes_KS        Same, but using a Kolmogorov-Smirnov-based 
#                   test. See docs for BiCopGofTest in the 
#                   VineCopula package
#GofRes_CvM_stat
#GofRes_KS_stat   Statistics for the tests
#Numboot          Number of bootstraps (numBSsmall/large)
#GofRes_Normal_CvM
#GofRes_Normal_KS
#GofRes_Normal_CvM_stat
#GofRes_Normal_KS_stat    Same as above, but for testing the
#                           goodness of fit of the normal 
#                           copula.
#Numboot_Normal   Number of bootstraps for the normal test 
#                   (numBSsmall/large)
#
OurBiCopSelect<-function(u1,u2,families,level=0.05,AICBIC="AIC",
                         numBSsmall=100,pthresh=0.2,numBSlarge=1000,
                         gofnormal=TRUE,status=TRUE)
{
  #first, test for independence (H0 is independence)
  if (status) {cat(paste("Starting independence test: ",Sys.time(),"\n"))}
  IndepTestRes<-BiCopIndTest(u1,u2)$p.value
  if (status) {cat(paste("Done:",Sys.time(),"\n"))}
  
  #if independence rejected, then get AICs for copulas and do 
  #goodness of fit stuff
  if (IndepTestRes<level){
    
    #AIC/BIC stuff
    if (status) {cat(paste("Starting A/BIC model selection: ",Sys.time(),"\n"))}
    InfCritRes<-data.frame(copcode=families,
                    copname=BiCopName(families, short=TRUE),
                    par1=NA,
                    par2=NA,
                    logLik=NA,
                    AIC=NA,
                    BIC=NA)
    for (counter in 1:(dim(InfCritRes)[1])){
      tres<-BiCopEst(u1,u2,family=InfCritRes[counter,1])
      InfCritRes$par1[counter]<-tres$par
      InfCritRes$par2[counter]<-tres$par2
      InfCritRes$logLik[counter]<-tres$logLik
      InfCritRes$AIC[counter]<-tres$AIC
      InfCritRes$BIC[counter]<-tres$BIC
    }
    if (status) {cat(paste("Done: ",Sys.time(),"\n"))}
    
    #g.o.f. stuff for the A/BIC-best copula
    if (status) {cat(paste("Starting gof for A/BIC-best copula: ",Sys.time(),"\n"))}
    if (AICBIC=="AIC"){
      ind<-which.min(InfCritRes$AIC)
    }
    if (AICBIC=="BIC"){
      ind<-which.min(InfCritRes$BIC)
    }
    if (AICBIC!="AIC" && AICBIC!="BIC"){
      stop("Error in OurBiCopSelect: incorrect AICBIC")
    }
    Numboot<-numBSsmall
    gres<-BiCopGofTest(u1,u2,family=InfCritRes$copcode[ind],
                       method="kendall",B=numBSsmall)
    GofRes_CvM<-gres$p.value.CvM
    GofRes_KS<-gres$p.value.KS
    if (GofRes_CvM<pthresh || GofRes_KS<pthresh)
    {
      Numboot<-numBSlarge
      gres<-BiCopGofTest(u1,u2,family=InfCritRes$copcode[ind],
                         method="kendall",B=numBSlarge)
      GofRes_CvM<-gres$p.value.CvM
      GofRes_KS<-gres$p.value.KS
    }
    GofRes_CvM_stat<-gres$statistic.CvM
    GofRes_KS_stat<-gres$statistic.KS
    if (status) {cat(paste("Done: ",Sys.time(),"\n"))}
    
    #g.o.f. stuff for the normal copula  
    if(gofnormal==T){
      if (status) {cat(paste("Starting gof for normal copula: ",Sys.time(),"\n"))}
      Numboot_Normal<-numBSsmall
      gres_normal_cop<-BiCopGofTest(u1,u2,family=1,method="kendall",B=numBSsmall)
      GofRes_Normal_CvM<-gres_normal_cop$p.value.CvM
      GofRes_Normal_KS<-gres_normal_cop$p.value.KS
      if (GofRes_Normal_CvM<pthresh || GofRes_Normal_KS<pthresh)
      {      
        Numboot_Normal<-numBSlarge
        gres_normal_cop<-BiCopGofTest(u1,u2,family=1,method="kendall",B=numBSlarge)
        GofRes_Normal_CvM<-gres_normal_cop$p.value.CvM
        GofRes_Normal_KS<-gres_normal_cop$p.value.KS
      }
      GofRes_Normal_CvM_stat<-gres_normal_cop$statistic.CvM
      GofRes_Normal_KS_stat<-gres_normal_cop$statistic.KS
      if (status) {cat(paste("Done: ",Sys.time(),"\n"))}
    }else{
      GofRes_Normal_CvM<-NA
      GofRes_Normal_KS<-NA
      GofRes_Normal_CvM_stat<-NA
      GofRes_Normal_KS_stat<-NA
      Numboot_Normal<-NA
    }
  } else {
    InfCritRes<-NA
    GofRes_CvM<-NA
    GofRes_KS<-NA
    GofRes_CvM_stat<-NA
    GofRes_KS_stat<-NA
    Numboot<-NA
    GofRes_Normal_CvM<-NA
    GofRes_Normal_KS<-NA
    GofRes_Normal_CvM_stat<-NA
    GofRes_Normal_KS_stat<-NA
    Numboot_Normal<-NA
  }
  
  return(list(IndepTestRes=IndepTestRes,
              InfCritRes=InfCritRes,
              GofRes_CvM=GofRes_CvM,
              GofRes_KS=GofRes_KS,
              GofRes_CvM_stat=GofRes_CvM_stat,
              GofRes_KS_stat=GofRes_KS_stat,
              Numboot=Numboot,
              GofRes_Normal_CvM=GofRes_Normal_CvM,
              GofRes_Normal_KS=GofRes_Normal_KS,
              GofRes_Normal_CvM_stat=GofRes_Normal_CvM_stat,
              GofRes_Normal_KS_stat=GofRes_Normal_KS_stat,
              Numboot_Normal=Numboot_Normal))
}


