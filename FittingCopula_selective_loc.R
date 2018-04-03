# THIS CODE TESTS FOR A GIVEN COPULA WHICH ARE THE POSSIBLE BEST FIT MODELS ?
#----------------
library(copula)
library(VineCopula)
#-------------------------------------
source("vivj_matrix.R")
source("good_loclist.R")
source(file="OurBiCopSelect.R")
#---------------------------------------------------------------------------------------
# This function generates result for single species
RES_single_sp<-function(sp,d_allsp,families,level,data_pt_thrs){
  
  len<-length(families)
  good_loc<-good_loclist(d_allsp=d_allsp,sp=sp,data_pt_thrs=data_pt_thrs)
  #cat("good_loc=",good_loc,"\n")
  
  lengoodloc<-length(good_loc)
  
  gfc_p_CvM<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_p_CvM) <- paste("loc",good_loc, sep="")
  rownames(gfc_p_CvM) <- paste("loc",good_loc, sep="")
  
  gfc_p_KS<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_p_KS)<-colnames(gfc_p_CvM)
  rownames(gfc_p_KS)<-rownames(gfc_p_CvM)
  
  gfc_p_CvM_stat<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_p_CvM_stat)<-colnames(gfc_p_CvM)
  rownames(gfc_p_CvM_stat)<-rownames(gfc_p_CvM)
  
  gfc_p_KS_stat<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_p_KS_stat)<-colnames(gfc_p_CvM)
  rownames(gfc_p_KS_stat)<-rownames(gfc_p_CvM)
  
  gfc_normal_p_CvM<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_normal_p_CvM)<-colnames(gfc_p_CvM)
  rownames(gfc_normal_p_CvM)<-rownames(gfc_p_CvM)
  
  gfc_normal_p_KS<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_normal_p_KS)<-colnames(gfc_p_CvM)
  rownames(gfc_normal_p_KS)<-rownames(gfc_p_CvM)
  
  gfc_normal_p_CvM_stat<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_normal_p_CvM_stat)<-colnames(gfc_p_CvM)
  rownames(gfc_normal_p_CvM_stat)<-rownames(gfc_p_CvM)
  
  gfc_normal_p_KS_stat<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_normal_p_KS_stat)<-colnames(gfc_p_CvM)
  rownames(gfc_normal_p_KS_stat)<-rownames(gfc_p_CvM)
  
  copdata_cor_Kend<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(copdata_cor_Kend)<-colnames(gfc_p_CvM)
  rownames(copdata_cor_Kend)<-rownames(gfc_p_CvM)
  
  gfc_numBS<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_numBS)<-colnames(gfc_p_CvM)
  rownames(gfc_numBS)<-rownames(gfc_p_CvM)
  
  gfc_numBS_success<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_numBS_success)<-colnames(gfc_p_CvM)
  rownames(gfc_numBS_success)<-rownames(gfc_p_CvM)
  
  gfc_normal_numBS<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_normal_numBS)<-colnames(gfc_p_CvM)
  rownames(gfc_normal_numBS)<-rownames(gfc_p_CvM)
  
  gfc_normal_numBS_success<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(gfc_normal_numBS_success)<-colnames(gfc_p_CvM)
  rownames(gfc_normal_numBS_success)<-rownames(gfc_p_CvM)
  
  LTdep_AICw<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(LTdep_AICw)<-colnames(gfc_p_CvM)
  rownames(LTdep_AICw)<-rownames(gfc_p_CvM)
  
  UTdep_AICw<-matrix(NA,nrow=lengoodloc,ncol=lengoodloc)
  colnames(UTdep_AICw)<-colnames(gfc_p_CvM)
  rownames(UTdep_AICw)<-rownames(gfc_p_CvM)
  
  info_ord_AIC<-array(NA,dim =c(lengoodloc,lengoodloc,len))  # To store all goodfit possibilities
  colnames(info_ord_AIC)<-colnames(gfc_p_CvM)
  rownames(info_ord_AIC)<-rownames(gfc_p_CvM)
  
  info_ord_copcode<-array(NA,dim =c(lengoodloc,lengoodloc,len))
  colnames(info_ord_copcode)<-colnames(gfc_p_CvM)
  rownames(info_ord_copcode)<-rownames(gfc_p_CvM)
  
  info_ord_copname<-array(NA,dim =c(lengoodloc,lengoodloc,len))
  colnames(info_ord_copname)<-colnames(gfc_p_CvM)
  rownames(info_ord_copname)<-rownames(gfc_p_CvM)
  
  info_ord_LTdep<-array(NA,dim =c(lengoodloc,lengoodloc,len)) 
  colnames(info_ord_LTdep)<-colnames(gfc_p_CvM)
  rownames(info_ord_LTdep)<-rownames(gfc_p_CvM)
  
  info_ord_UTdep<-array(NA,dim =c(lengoodloc,lengoodloc,len))
  colnames(info_ord_UTdep)<-colnames(gfc_p_CvM)
  rownames(info_ord_UTdep)<-rownames(gfc_p_CvM)
  
  info_ord_AICw<-array(NA,dim =c(lengoodloc,lengoodloc,len))
  colnames(info_ord_AICw)<-colnames(gfc_p_CvM)
  rownames(info_ord_AICw)<-rownames(gfc_p_CvM)
  
 
  
#  cat("sp=",sp,"\n")
  
  num_indep<-0
  num_neg_cor<-0
  
  for(i in c(1:lengoodloc)){
    for(j in c(1:lengoodloc)){
      if(i!=j){
        m<-vivj_matrix(d_allsp,sp,good_loc[i],good_loc[j])
        u1<-m[,1]
        u2<-m[,2]
        #plot(u1,u2)
#        cat("[i, j]=[",i," ,",j,"]","(good_loc[i],good_loc[j])=(",good_loc[i],",",good_loc[j],")","\n")
        
        ans<-OurBiCopSelect(u1,u2,families,level=0.05,AICBIC="AIC",
                            numBSsmall=100,pthresh=0.2,numBSlarge=1000,
                            gofnormal=FALSE,status=FALSE)
        
        copdata_cor_Kend[i,j]<-ans$TauVal
        
        if(ans$IndepTestRes<level && ans$TauVal>0){
          gfc_numBS[i,j]<-ans$Numboot
          gfc_numBS_success[i,j]<-ans$Numboot_success
          gfc_p_CvM[i,j]<-ans$GofRes_CvM
          gfc_p_KS[i,j]<-ans$GofRes_KS
          gfc_p_CvM_stat[i,j]<-ans$GofRes_CvM_stat
          gfc_p_KS_stat[i,j]<-ans$GofRes_KS_stat
          
          gfc_normal_numBS[i,j]<-ans$Numboot_Normal
          gfc_normal_numBS_success[i,j]<-ans$Numboot_success_Normal
          gfc_normal_p_CvM[i,j]<-ans$GofRes_Normal_CvM
          gfc_normal_p_KS[i,j]<-ans$GofRes_Normal_KS
          gfc_normal_p_CvM_stat[i,j]<-ans$GofRes_Normal_CvM_stat
          gfc_normal_p_KS_stat[i,j]<-ans$GofRes_Normal_KS_stat
          
          LTdep_AICw[i,j]<-ans$relLTdep_AICw
          UTdep_AICw[i,j]<-ans$relUTdep_AICw
          
          info<-ans$InfCritRes
          info_ord<-info[order(info[,6],decreasing = F),] # ordered according to min to max AIC
          info_ord_AIC[i,j,]<-info_ord$AIC
          info_ord_copcode[i,j,]<-info_ord$copcode
          info_ord_copname[i,j,]<-as.character(info_ord$copname)
          info_ord_LTdep[i,j,]<-info_ord$LTdep
          info_ord_UTdep[i,j,]<-info_ord$UTdep
          info_ord_AICw[i,j,]<-info_ord$AICw
 
        }else if(ans$IndepTestRes<level && ans$TauVal<0){
          num_neg_cor<-num_neg_cor+1
          gfc_numBS[i,j]<- -Inf  # I just put -Inf to see when it's -vely correlated data?
          gfc_numBS_success[i,j]<- -Inf
          gfc_p_CvM[i,j]<- -Inf         
          gfc_p_KS[i,j]<- -Inf
          gfc_p_CvM_stat[i,j]<- -Inf         
          gfc_p_KS_stat[i,j]<- -Inf
          gfc_normal_numBS[i,j]<- -Inf
          gfc_normal_numBS_success[i,j]<- -Inf
          gfc_normal_p_CvM[i,j]<- -Inf  
          gfc_normal_p_KS[i,j]<- -Inf
          gfc_normal_p_CvM_stat[i,j]<- -Inf
          gfc_normal_p_KS_stat[i,j]<- -Inf
          info_ord_AIC[i,j,]<- -Inf
          info_ord_copcode[i,j,]<- -Inf
          info_ord_copname[i,j,]<-"NegCor"
          info_ord_LTdep[i,j,]<- -Inf
          info_ord_UTdep[i,j,]<- -Inf
          info_ord_AICw[i,j,]<- -Inf
          LTdep_AICw[i,j]<- -Inf
          UTdep_AICw[i,j]<- -Inf
          
        }else{
          num_indep<-num_indep+1
          gfc_numBS[i,j]<-Inf
          gfc_numBS_success[i,j]<-Inf
          gfc_p_CvM[i,j]<-Inf         # I just put Inf to see when it's indep?
          gfc_p_KS[i,j]<-Inf
          gfc_p_CvM_stat[i,j]<-Inf         # I just put Inf to see when it's indep?
          gfc_p_KS_stat[i,j]<-Inf
          gfc_normal_numBS[i,j]<-Inf
          gfc_normal_numBS_success[i,j]<-Inf
          gfc_normal_p_CvM[i,j]<-Inf  # If gofnormal==F, then gfc_rmal_p matrices also contains Inf at all off-diagonal entries
          gfc_normal_p_KS[i,j]<-Inf
          gfc_normal_p_CvM_stat[i,j]<-Inf
          gfc_normal_p_KS_stat[i,j]<-Inf
          info_ord_AIC[i,j,]<-Inf
          info_ord_copcode[i,j,]<-Inf
          info_ord_copname[i,j,]<-"Indep"
          info_ord_LTdep[i,j,]<-Inf
          info_ord_UTdep[i,j,]<-Inf
          info_ord_AICw[i,j,]<-Inf
          LTdep_AICw[i,j]<-Inf
          UTdep_AICw[i,j]<-Inf
        }
        
      }
    }
    
  }
  #-----------------------------------------------------------------------------------------------
  # Save the results
  RES_sp<-list(num_indep_loc_pair=num_indep/2,
               num_neg_cor_loc_pair=num_neg_cor/2,
               gfc_numBS=gfc_numBS,
               gfc_numBS_success=gfc_numBS_success,
               gfc_p_CvM=gfc_p_CvM,
               gfc_p_KS=gfc_p_KS,
               gfc_p_CvM_stat=gfc_p_CvM_stat,
               gfc_p_KS_stat=gfc_p_KS_stat,
               gfc_normal_numBS=gfc_normal_numBS,
               gfc_normal_numBS_success=gfc_normal_numBS_success,
               gfc_normal_p_CvM=gfc_normal_p_CvM,
               gfc_normal_p_KS=gfc_normal_p_KS,
               gfc_normal_p_CvM_stat=gfc_normal_p_CvM_stat,
               gfc_normal_p_KS_stat=gfc_normal_p_KS_stat,
               info_ord_copcode=info_ord_copcode,
               info_ord_copname=info_ord_copname,
               info_ord_AIC=info_ord_AIC,
               info_ord_LTdep=info_ord_LTdep,
               info_ord_UTdep=info_ord_UTdep,
               info_ord_AICw=info_ord_AICw,
               LTdep_AICw=LTdep_AICw,
               UTdep_AICw=UTdep_AICw,
               copdata_cor_Kend=copdata_cor_Kend)
  
  return(RES_sp)
  
} 

#--------------------------------------------------------------------------











