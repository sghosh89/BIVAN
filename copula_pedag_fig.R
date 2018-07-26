library(VineCopula)
# This function takes 
# Input to create a seed copula:
#        fam : family of considered bivariate copula with two parameters
#        p1range : range of par
#        p2range : range of par2
#        incr : increment for parameter range
# and gives 
# Output after BiCopEst of that specified copula:
#       a data.frame which gives info like estimated par,par2,kendall's tau,taildep

cop_est_2par<-function(fam,p1range,p2range,incr){
  par_range<-as.vector(seq(from=p1range[1],to=tail(p1range,1),by=incr))
  par2_range<-as.vector(seq(from=p2range[1],to=tail(p2range,1),by=incr))
  dm<-length(par_range)*length(par2_range)
  df<-data.frame(par=NA*numeric(dm),par2=NA*numeric(dm),
                 Tau=NA*numeric(dm),
                 LT=NA*numeric(dm),
                 UT=NA*numeric(dm),
                 LTmUT=NA*numeric(dm))
  count<-0
  for(par in par_range){
    for(par2 in par2_range){
      count<-count+1
      #cat("par=",par,"par2=",par2,"\n")
      cop<-BiCop(family=fam,par=par,par2=par2)
      df$par[count]<-cop$par
      df$par2[count]<-cop$par2
      df$Tau[count]<-cop$tau
      df$LT[count]<-cop$taildep$lower 
      df$UT[count]<-cop$taildep$upper 
      df$LTmUT[count]<-(cop$taildep$lower - cop$taildep$upper)
    }
  }
  
  return(df)
}
#-----------------------------------------------
# This function takes input :
#            df : A dataframe which is the output of cop_est_2par function
#            ncut : an integer which indicates the number of label in LTmUT range
#            Tau0 : target Tau with which you want to plot your copula
#            eps : error limit in Tau0
#            LTmUT_vector_given : logical (if F then it will generate LTmUT_vector automatically from df)
#            my_LTmUT_vector : a vector for LTmUT is given if the above argument is T, otherwise it's NA

# gives output is a list of two:
#            A screened version of input dataframe with specific LTmUT dep value
#            LTmUT_vector
screening_df<-function(df,ncut,Tau0,eps,LTmUT_vector_given,my_LTmUT_vector){
  
  ind<-which(abs(df$Tau-Tau0)<eps)
  df<-df[ind,]
  
  if(LTmUT_vector_given==F){
    LmU<-range(df$LTmUT)
    LTmUT<-trunc(min(abs(LmU))*10)/10 # truncating upto 1st decimal place
    ndiv<-ncut-1
    LTmUT_vector<-seq(from=-LTmUT,to=LTmUT,by=(2*LTmUT)/ndiv)
  }else{
    LTmUT_vector<-my_LTmUT_vector
  }
  
  id<-c()
  for(i in c(1:length(LTmUT_vector))){
    id<-c(id,which.min(abs(df$LTmUT - LTmUT_vector[i])))
  }
  df_screened<-df[id,]
  return(list(df_screened=df_screened,
              LTmUT_vector=LTmUT_vector))
}
#----------------------

