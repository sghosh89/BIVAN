#A function for getting information on Taylor's law, both spatial and temporal
#
#Args
#m    A matrix, time series down the columns. Assumed no 0s or NAs, no variances 
#       are zero computed either down the columns or across the rows
#
#Output - a list containing:
#1) an element named s.stats which is a named vector with these named elements
#s.lin      p-value result for test of linearity of spatial TL
#s.hom      p-value result for test of homoskedasticity of spatial TL
#s.rmse     the root mean squared error of the spatial TL model
#s.int      intercept of the same relationship
#s.slope    slope of the spatial log(v)-log(m) relationships
#2) an element named t.stats which is a named vector with these named elements
#t.lin      this and the following similar to the above but for temporal TL
#t.hom
#t.rmse
#t.int      
#t.slope
#3) an element named s.points which is the log10 means and log10 variances
#for spatial TL
#4) similar for temporal
#
tl_stats<-function(m)
{
  sres<-worker_tl(m)
  s.stats<-sres$stats
  names(s.stats)<-c("s.lin","s.hom","s.rmse","s.int","s.slope","s.quad.coeff")
  s.points<-sres$points
  
  tres<-worker_tl(t(m))
  t.stats<-tres$stats
  names(t.stats)<-c("t.lin","t.hom","t.rmse","t.int","t.slope","t.quad.coeff")
  t.points<-tres$points
  
  return(list(s.stats=s.stats,t.stats=t.stats,
              s.points=s.points,t.points=t.points))
}

worker_tl<-function(m)
{
  N<-dim(m)[1]
  n<-dim(m)[2]
  
  lms<-log10(apply(X=m,MARGIN=1,FUN=mean))
  lvs<-log10(apply(X=m,MARGIN=1,FUN=var))
  
  lms2<-lms^2
  mod1<-lm(lvs~lms)
  mod2<-lm(lvs~lms2+lms)
  hlin<-anova(mod1,mod2) #make the quadratic comparison
  absresids<-abs(resid(mod1))
  predicts<-predict(mod1)
  mod3<-lm(absresids~predicts)
  hhet<-anova(mod3)
  
  res<-c(unname(hlin[6][2,1]), #the p-value for the quadratic comparison
         unname(hhet[5][1,1]), #the p-value for a test of homoskedasticity
         sqrt(mean(resid(mod1)^2)), #the root mean squared error for the TL model
         unname(mod1$coefficients[1]), #the intercept
         unname(mod1$coefficients[2]), #the slope
         unname(mod2$coefficients[2]))  # the coefficient of quadratic term
  
  return(list(stats=res,points=data.frame(log10m=lms,log10v=lvs)))
}