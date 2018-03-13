


#This function takes a bunch of time series measured at the same
#times in different locations and creates surrogate datasets ...
#
#This differs from copsurrog2d and ncsurrog because...
#
#Args
#m          A N by n matrix, where N is the length of the time series 
#           and n is the number of time series
#targetcop  The target copula. This is a copula object from the 
#             copula package. Must have dimension equal to dim(m)[2]

#Output
#

copsurrognd<-function(m,tcop)
{
  N<-dim(m)[1]
  n<-dim(m)[2]
  
  #error checking
  #check the dimension of targetcop is right

  
}