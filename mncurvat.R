#Gives the average curvature of a quadratic function f(x)=ax^2+bx+c
#at a bunch of values of x
#
#coefs      c(a,b,c)
#xvals      The values of x
#
#Returns one single number which is the average curvature

mncurvat<-function(coefs,xvals){
  a<-coefs[1]
  b<-coefs[2]
  c<-coefs[3]
  kappa<-2*a/((1+(2*a*xvals+b)^2)^(3/2))
  return(mean(kappa))
}

