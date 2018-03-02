#-----------------------------------------------------------------------------------------------------------------------
# THIS CODE CONTAINS ALL THE STATISTISCAL FUNCTIONS WHICH CAN BE TESTED ON THE TWO RANK VECTORS OF A GIVEN COPULA
#-----------------------------------------------------------------------------------------------------------------------

#--------------------------- STATISTICS 2 : correlation based Stat ---------------------------------------------------------

source("CopulaFunction_flexible.R")

CorlCoru<-function(vi,vj)
{
  Corl<-Corbds(vi,vj,0,0.5)
  Coru<-Corbds(vi,vj,0.5,1)

  return(c(Corl,Coru))
}


#--------------------------- STATISTICS 4 : Stat based on the distance from right diagonal in lower and upper triangle ---------------------------------------------------------

PlPu<-function(vi,vj){
  
  ind_ll<-which(vi+vj<1)
  ind_ur<-which(vi+vj>1)
  
  if(length(ind_ll)!=0 & length(ind_ur)!=0){ 
    
    distance_ll<-sort(abs(vi[ind_ll]-vj[ind_ll])/sqrt(2))
    distance_ur<-sort(abs(vi[ind_ur]-vj[ind_ur])/sqrt(2))
    d_max<-sqrt(2)/2
    
    Sl_P_c<-0
    Sl_P_uniq<-c(0)
    distance_ll_df<-as.data.frame(table(distance_ll))
    for (i in 1:length(distance_ll_df$Freq)){
      Sl_P_c<-Sl_P_c+(distance_ll_df$Freq[i]/length(distance_ll))
      Sl_P_uniq<-c(Sl_P_uniq,Sl_P_c)
    }
    Sl_P_uniq<-c(Sl_P_uniq,1)
    dist_Sl_P_uniq<-c(0,as.numeric(as.character(distance_ll_df$distance_ll)),d_max)
    
    Su_P_c<-0
    Su_P_uniq<-c(0)
    distance_ur_df<-as.data.frame(table(distance_ur))
    for (i in 1:length(distance_ur_df$Freq)){
      Su_P_c<-Su_P_c+(distance_ur_df$Freq[i]/length(distance_ur))
      Su_P_uniq<-c(Su_P_uniq,Su_P_c)
    }
    Su_P_uniq<-c(Su_P_uniq,1)
    dist_Su_P_uniq<-c(0,as.numeric(as.character(distance_ur_df$distance_ur)),d_max)
    
    nrep<-2
    Sl_P<-rep(Sl_P_uniq,each=nrep)
    dist_Sl_P<-rep(dist_Sl_P_uniq[2:length(dist_Sl_P_uniq)],each=nrep)
    dist_Sl_P<-c(0,dist_Sl_P,d_max)
    
    Su_P<-rep(Su_P_uniq,each=nrep)
    dist_Su_P<-rep(dist_Su_P_uniq[2:length(dist_Su_P_uniq)],each=nrep)
    dist_Su_P<-c(0,dist_Su_P,d_max)
    
    dist_Si_P<-c()
    Si_P<-c()
    for (di in seq(from=0, to=d_max, by=d_max/1000)){
      dist_Si_P<-c(dist_Si_P,di)
      a1<-(1/sqrt(2))- di
      Si_P<-c(Si_P,2*di*(a1+(1/sqrt(2))))
    } 
    
    Sl_Su_Si_P<-list(dist_Sl_P=dist_Sl_P,Sl_P=Sl_P,dist_Su_P=dist_Su_P,Su_P=Su_P,dist_Si_P=dist_Si_P,Si_P=Si_P)
    
    #integrals under the functions Sl_P, Su_P, Si_P
    Au_Sl_P<-0
    Au_Su_P<-0 
    
    for(i in 1:(length(Sl_P_uniq)-1)){
      Au_Sl_P<-Au_Sl_P+(Sl_P_uniq[i]*(dist_Sl_P_uniq[i+1]-dist_Sl_P_uniq[i]))
    }
    
    for(i in 1:(length(Su_P_uniq)-1)){
      Au_Su_P<-Au_Su_P+(Su_P_uniq[i]*(dist_Su_P_uniq[i+1]-dist_Su_P_uniq[i]))
    }
    
    Au_Si_P<-((d_max^2)*sqrt(2))-(2*(d_max^3)/3)   # Analytical integration
    
    Pl<-Au_Sl_P-Au_Si_P
    Pu<-Au_Su_P-Au_Si_P
    
  }else{
    Sl_Su_Si_P<-list(dist_Sl_P=NA,Sl_P=NA,dist_Su_P=NA,Su_P=NA,dist_Si_P=NA,Si_P=NA)
    Pl<-NA
    Pu<-NA
  }
  
  return(list(Sl_Su_Si_P,Pl,Pu)) 
}

#---------------- STATISTICS 5 : stat based on the distance for lower and upper triangle of left diagonal of the box -----------------------------------------------------
# d_max is the distance from corner to the diagonal of the box

TlTu<-function(vi,vj){
  
  i_lt<-which(vi+vj<1)
  i_ut<-which(vi+vj>1)
  
  if(length(i_lt)!=0 & length(i_ut)!=0){ 
    
    distance_lt<-sort(abs(vi[i_lt]+vj[i_lt])/sqrt(2))
    distance_ut<-sort(abs(vi[i_ut]+vj[i_ut]-2)/sqrt(2))
    d_max<-sqrt(2)/2
    
    Il_c<-0
    Il_uniq<-c(0)
    distance_lt_df<-as.data.frame(table(distance_lt))
    for (i in 1:length(distance_lt_df$Freq)){
      Il_c<-Il_c+(distance_lt_df$Freq[i]/length(distance_lt))
      Il_uniq<-c(Il_uniq,Il_c)
    }
    Il_uniq<-c(Il_uniq,1)
    dist_Il_uniq<-c(0,as.numeric(as.character(distance_lt_df$distance_lt)),d_max)
    
    Iu_c<-0
    Iu_uniq<-c(0)
    distance_ut_df<-as.data.frame(table(distance_ut))
    for (i in 1:length(distance_ut_df$Freq)){
      Iu_c<-Iu_c+(distance_ut_df$Freq[i]/length(distance_ut))
      Iu_uniq<-c(Iu_uniq,Iu_c)
    }
    Iu_uniq<-c(Iu_uniq,1)
    dist_Iu_uniq<-c(0,as.numeric(as.character(distance_ut_df$distance_ut)),d_max)
    
    nrep<-2
    Il<-rep(Il_uniq,each=nrep)
    dist_Il<-rep(dist_Il_uniq[2:length(dist_Il_uniq)],each=nrep)
    dist_Il<-c(0,dist_Il,d_max)
    
    Iu<-rep(Iu_uniq,each=nrep)
    dist_Iu<-rep(dist_Iu_uniq[2:length(dist_Iu_uniq)],each=nrep)
    dist_Iu<-c(0,dist_Iu,d_max)
    
    dist_Ii<-c()
    Ii<-c()
    numpts<-1000
    for (d_th in seq(from=0, to=d_max, by=d_max/numpts)){
      dist_Ii<-c(dist_Ii,d_th)
      Ii<-c(Ii,(sqrt(2)*d_th)^2)
    } 
    
    Il_Iu_Ii<-list(dist_Il=dist_Il,Il=Il,dist_Iu=dist_Iu,Iu=Iu,dist_Ii=dist_Ii,Ii=Ii)
    
    #integrals under the functions Il, Iu, Ii
    Au_Il<-0
    Au_Iu<-0 
    
    for(i in 1:(length(Il_uniq)-1)){
      Au_Il<-Au_Il+(Il_uniq[i]*(dist_Il_uniq[i+1]-dist_Il_uniq[i]))
    }
    
    for(i in 1:(length(Iu_uniq)-1)){
      Au_Iu<-Au_Iu+(Iu_uniq[i]*(dist_Iu_uniq[i+1]-dist_Iu_uniq[i]))
    }
    
    Au_Ii<-(2/3)*(d_max^3)   # Analytical integration
    
    Tl<-Au_Il-Au_Ii
    Tu<-Au_Iu-Au_Ii
    
  }else{
    Il_Iu_Ii<-list(dist_Il=NA,Il=NA,dist_Iu=NA,Iu=NA,dist_Ii=NA,Ii=NA)
    Tl<-NA
    Tu<-NA
  }
  
  return(list(Il_Iu_Ii,Tl,Tu))
  
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