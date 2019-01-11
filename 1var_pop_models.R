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
Ricker(r=1.1,K=50,p0=5,lensim=100) # bifurcation point r=1
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
  
  #provided b>2 (to get finite r2), these are the stability conditions:
  cat("For stable, monotonic approach towards eqm. upper limit of r should be",r1,"\n")

  cat("For stable, oscillatory approach towards eqm. upper limit of r should be within ",r1,"to",r2,"\n")
  
}

#call the function
r=5
a=0.5
b=1.5
K_e<-((r^(1/b))-1)/a
K_e
Hassell(r=r, a=a, b=b, p0=K_e, lensim=200) #b~1.5 : contest, b~100 : scramble

Hassell(r=5, a=0.01, b=2.5, p0=0.1, lensim=200) # see change : r=5,a=0.5 and vary b 1.5 to 5.5
#----------------------------------------------------------------------------------------------------
# Exploring Maynard Smith model
Msmith<-function(r,a,b,p0,lensim){
  pop<-c(p0)
  time<-c(0)
  
  for(it in c(1:lensim)){
    pt<-(r*p0)/(1+((a*p0)^b))
    time<-c(time,it)
    pop<-c(pop,pt)
    p0<-pt
  }
  
  plot(time,pop,type="b")
  
}

#call the function
Msmith(r=4.5,a=0.5,b=1.4,p0=0.1,lensim=200) #see changes as r=5,a=0.5 and vary b : 1,1.5,2,4
                                          #see changes as a=0.5,b=4 and vary r : 1.2, 1.5, 2





