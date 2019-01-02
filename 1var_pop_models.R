# Exploring Ricker model
Ricker<-function(r,K,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  # po<-rep(K,r)
  # print(po)
  
  for(it in c(1:lensim)){
    pt<-p0*exp(r*(1-(p0/K)))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
Ricker(r=2,K=100,p0=25,lensim=100)
#-------------------------------------------------------

# Exploring Hassell model
# help: https://jmahaffy.sdsu.edu/courses/s00a/math121/lectures/qual_discrete/qualdiscrete.html#Hassell
#       https://www.jstor.org/stable/pdf/3863.pdf?refreqid=excelsior%3A348907abb1bca7f039109d94db5667b5
#       http://lab.rockefeller.edu/cohenje/PDFs/231CohenUnexpectedDominanceHighFreqChaoticNonlinPpnModelsLet.pdf

Hassell<-function(r,a,b,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-(r*p0)/((1+(a*p0))^b)
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  r1<-((b/(b-1)))^b
  r2<-((b/(b-2)))^b
  
  plot(time,pop,type="b")
  cat("For stable, monotonic approach towards eqm. upper limit of r should be",r1,"\n")
  
  cat("For stable, oscillatory approach towards eqm. lower limit of r should be",r1,"\n")
  cat("For stable, oscillatory approach towards eqm. upper limit of r should be",r2,"\n")
  
}

#call the function
Hassell(r=5, a=0.01, b=4, p0=25, lensim=200)

#----------------------------------------------------------------------------------------------------
# Exploring Maynard Smith model







